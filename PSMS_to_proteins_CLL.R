library(igraph)
library(tidyr)
library(dplyr)
library(vsn)
library(parallel)
'%!in%' <- function(x,y)!('%in%' (x,y))
expand_rows_with_multiple_terms <- function(df, column_name) {
  
  df <- df %>%
    mutate(!!column_name := strsplit(as.character(.[[column_name]]), ",")) %>%
    unnest(!!column_name)
  
  return(df)
}
randomize_data <- function(data) {
  if("membership" %in% names(data)) {
    #shuffle only the membership column
    data$membership <- sample(data$membership)
  }
  return(data)
}


setwd("~/Documents/githubspot/ibrutinib-cll/")

cll_psms <- read.delim("~/Documents/githubspot/ibrutinib-cll/data/cll_psms.txt") 

#center and scale
normalized_sets <- lapply(split(cll_psms, cll_psms$Biological.set), function(set_data) {
  tmt_data <- as.matrix(select(set_data, starts_with("tmt10plex")))
  vsn_normalized <- exprs(vsn2(tmt_data, returnData = TRUE))
  bind_cols(select(set_data, -starts_with("tmt10plex")), as.data.frame(vsn_normalized))
})
cll_psms <- bind_rows(normalized_sets)

graphs_comms <- readRDS("~/Documents/githubspot/ibrutinib-cll/data/graphs_comms_again.RDS")

hits_everything <- readRDS("~/Documents/githubspot/ibrutinib-cll/output/hits_everything.RDS")



cll_meta <- read.delim("~/Downloads/CLL_metadata 1.txt") 
cll_meta <- expand_rows_with_multiple_terms(cll_meta, "TMT")

get_CLL_and_annotate <- function(CLL_data_psms, graphs_comms, ioi){peptide_seq_list <- list()


peptide_seq_list[[ioi]] <- data.frame(
  
  as.data.frame(vertex_attr(graphs_comms[[ioi]])[["name"]]),
  # as.data.frame(gsub("[0-9\\.\\+]", "",
  #                    vertex_attr(graphs_comms[[i]])[["name"]])),
  
  as.data.frame(vertex_attr(graphs_comms[[ioi]])[["membership"]]))
colnames(peptide_seq_list[[ioi]]) <- c("peptides", #"peptides_full", 
                                       "membership" )


transform_data <- function(data) {
  data <- gsub("(\\w)\\+15\\.995", "\\L\\1", data, perl = TRUE)
  data <- gsub("(\\w)\\+57\\.021", "\\L\\1", data, perl = TRUE)
  data <- gsub("\\+\\d+\\.\\d+", "", data)
  
  return(data)
}

#transform peptides
peptide_seq_list[[ioi]]$peptides <- transform_data(data = peptide_seq_list[[ioi]]$peptides)
CLL_data_psms$Peptide <- transform_data(data = CLL_data_psms$Peptide)
#filter PSM file by ioi
CLL_data_psms_filt <- CLL_data_psms[CLL_data_psms$Peptide %in% peptide_seq_list[[ioi]]$peptides,]

CLL_data_psms_filt$membership <- 0

for(i in CLL_data_psms_filt$Peptide){
  
  CLL_data_psms_filt[CLL_data_psms_filt$Peptide == i,]$membership <- 
    peptide_seq_list[[ioi]][peptide_seq_list[[ioi]]$peptides == i,]$membership
  
}

#label sample source
CLL_data_psms_filt <- CLL_data_psms_filt %>%
  gather(key = "tmt_key", value = "tmt_value", starts_with("tmt10plex")) %>%
  mutate(type = paste(Biological.set, tmt_key, sep = "_"))%>%
  select(-tmt_key, everything())



return(CLL_data_psms_filt)

}

#test function 
CLL_data_psms_filt <- get_CLL_and_annotate(cll_psms, graphs_comms, ioi = "BTK")

#Create and process items to prep for the new data
CLL_only_hngc <- unique(cll_psms$Gene.Symbol[cll_psms$Gene.Symbol %!in% hits_everything$hgnc_symbol &
                                               !grepl(";", cll_psms$Gene.Symbol)])
all_from_CLL <- data.frame("id" = hits_everything$id, "gene symbol" = hits_everything$hgnc_symbol,
                           "test" = hits_everything$type, "pAdj" = hits_everything$pAdj)

CLL_ids_hngc <- data.frame("id" = "N.D.", "gene symbol" = CLL_only_hngc,
                           "test" = "CLL only", "pAdj" = "N.D.")
all_from_CLL <- rbind(all_from_CLL, CLL_ids_hngc)

merge_columns <- c(outer(unique(cll_psms$Biological.set), colnames(cll_psms)[grepl("tmt", colnames(cll_psms))], FUN = function(x, y) paste0(x, "_", y)))


for (new_col_name in merge_columns) {
  all_from_CLL[[new_col_name]] <- "N.D"
}

all_from_CLL <-all_from_CLL[!is.na(all_from_CLL$id),]

#parallel processing function
process_subset <- function(gene_symbol, membership, cll_psms, graphs_comms, all_from_CLL) {
  subset_id <- paste0(gene_symbol, "_", membership)
  
  #collect data once
  CLL_data_psms_filt <- tryCatch(
    get_CLL_and_annotate(CLL_data_psms = cll_psms, graphs_comms = graphs_comms, ioi = gene_symbol),
    error = function(e) { data.frame() }  #empty data frame on error
  )
  
  #store modifications for each mode
  modifications_list <- list()
  
  #process data based on mode
  process_data <- function(data, mode) {
    modifications <- list()
    if (nrow(data) > 0) {
      if (mode == "gene_symbol") {
        data$membership <- 0 #membership = 0 for gene symbol mode
      }
      
      all_types <- unique(data$type)
      for (type in all_types) {
        type_specific_data <- data[data$type == type , ]
        if(mode != "gene_symbol"){
          
          type_specific_data <- type_specific_data[type_specific_data$membership == membership , ]
        }
        if("tmt_value" %in% names(type_specific_data) && !is.null(type_specific_data$tmt_value) && length(type_specific_data$tmt_value) > 0) {
          type_sum <- sum(type_specific_data$tmt_value, na.rm = TRUE)
          
          if (type_sum != 0) {
            modifications[[type]] <- type_sum
          }
        }
      }
    }
    return(modifications)
  }
  
  #process default mode
  modifications_list[["default"]] <- process_data(CLL_data_psms_filt, "default")
  
  #process gene symbol mode
  modifications_list[["gene_symbol"]] <- process_data(CLL_data_psms_filt, "gene_symbol")
  
  #process random mode
  #randomize data first
  randomized_data <- randomize_data(CLL_data_psms_filt)
  modifications_list[["random"]] <- process_data(randomized_data, "random")
  
  #Return: a list with id, modifications for each mode
  return(list(id = subset_id, modifications = modifications_list))
}



#apply the function
unique_ids <- unique(all_from_CLL$id)

# for (id in unique_ids) {
#   #split into gene_symbol and membership
#   parts <- strsplit(id, "_", fixed = TRUE)[[1]]
#   gene_symbol <- parts[1]
#   membership <- parts[2]
#   
#   all_from_CLL <- process_subset(gene_symbol, membership, cll_psms, graphs_comms, all_from_CLL)
# }

#set up parallel cores
numCores <- detectCores() - 1
cl <- makeCluster(numCores)

#Run in parallel cluster
clusterExport(cl, list("process_subset", "cll_psms","randomize_data", "get_CLL_and_annotate", "graphs_comms", "all_from_CLL", "unique_ids"))
clusterEvalQ(cl, {
  library(dplyr)
  library(tidyr)
  library(igraph)
})

results <- parLapply(cl, unique_ids, function(id) {
  parts <- strsplit(id, "_", fixed = TRUE)[[1]]
  gene_symbol <- parts[1]
  membership <- parts[2]
  
  process_subset(gene_symbol, membership, cll_psms, graphs_comms, all_from_CLL)
})

columns_needed <- names(all_from_CLL)
all_from_CLL_default <- all_from_CLL[0, ]  #empty dataframe with the same columns
all_from_CLL_gene_symbol <- all_from_CLL[0, ]
all_from_CLL_random <- all_from_CLL[0, ]


  
  #Update a specific mode dataframe
  update_mode_df <- function(mode_df, mode_modifications, id, all_from_CLL) {
    #Find the row indices that match the id
    row_indices <- which(mode_df$id == id)
    
    #If the ID is not found, add a new row
    if (length(row_indices) == 0) {
      #Find the row that matches the id
      source_row_index <- which(all_from_CLL$id == id)
      
      if (length(source_row_index) > 0) {
        #Add the row to the mode-specific dataframe
        new_row <- all_from_CLL[source_row_index, , drop = FALSE]
        mode_df <- rbind(mode_df, new_row)
        row_indices <- (nrow(mode_df) - nrow(new_row) + 1):nrow(mode_df)
      } else {
        #If not found, create a new row with the same columns
        new_row <- setNames(as.data.frame(matrix(ncol = length(names(all_from_CLL)), nrow = 1)), names(all_from_CLL))
        new_row$id <- id
        mode_df <- rbind(mode_df, new_row)
        row_indices <- nrow(mode_df)
      }
    }
    
    #Update all matching rows
    for (col_name in names(mode_modifications)) {
      if(col_name %in% names(mode_df)) {
        mode_df[row_indices, col_name] <- mode_modifications[[col_name]]
      }
    }
    
    return(mode_df)
  }
  

for (result in results) {
  id <- result$id
  modifications_list <- result$modifications
  
  #Call the update function for each mode without changing the id
  all_from_CLL_default <- update_mode_df(all_from_CLL_default, modifications_list[["default"]], id, all_from_CLL)
  all_from_CLL_gene_symbol <- update_mode_df(all_from_CLL_gene_symbol, modifications_list[["gene_symbol"]], id, all_from_CLL)
  all_from_CLL_random <- update_mode_df(all_from_CLL_random, modifications_list[["random"]], id, all_from_CLL)
}


#Stop the cluster
stopCluster(cl)
###############

#Normalize and center
normalize_and_center <- function(df) {
  #Once numeric, then log2-transform
  df_normcenter <- df %>%
    mutate(across(starts_with("Set"), ~ as.numeric(.))) %>%
    mutate(across(starts_with("Set"), ~ log2(. + 1))) # Adding 1 to avoid log2(0)
  
  #Center by subtracting the median of each row
  sample_columns <- grep("^Set", names(df_normcenter), value = TRUE)
  df_normcenter[sample_columns] <- sweep(df_normcenter[sample_columns], 1, apply(df_normcenter[sample_columns], 1, median, na.rm = TRUE), "-")
  
  #Center by subtracting the column median
  df_normcenter <- df_normcenter %>%
    mutate(across(starts_with("Set"), ~ . - median(., na.rm = TRUE)))
  
  return(df_normcenter)
}


all_from_CLL_default_normcenter <- normalize_and_center(all_from_CLL_default)
all_from_CLL_gene_symbol_normcenter <- normalize_and_center(all_from_CLL_gene_symbol)
all_from_CLL_random_normcenter <- normalize_and_center(all_from_CLL_random)


###############
#Save
first_all_from_CLL_normcenter <- readRDS("data/re_search_all_from_CLL.RDS")
saveRDS(all_from_CLL_default_normcenter, "data/all_from_CLL_default_normcenter.RDS")
saveRDS(all_from_CLL_random_normcenter, "data/all_from_CLL_random_normcenter.RDS")
saveRDS(all_from_CLL_gene_symbol_normcenter, "data/all_from_CLL_gene_symbol_normcenter.RDS")

###############
#Randomize again per peptide
randomize_data2 <- function(data) {
  if("Peptide" %in% names(data)){
    data <- as_tibble(data)
    
    peptide_groups <- unique(data$Peptide)
    
    memberships <- data %>%
      distinct(Peptide, membership) %>%
      pull(membership)
    
    shuffled_memberships <- sample(memberships)
    
    lookup_table <- tibble(Peptide = peptide_groups, membership = shuffled_memberships)
    
    data <- data %>%
      select(-membership) %>%
      left_join(lookup_table, by = "Peptide")
  }
  
  return(as.data.frame(data))
}

process_subset2 <- function(gene_symbol, membership, cll_psms, graphs_comms, all_from_CLL) {
  subset_id <- paste0(gene_symbol, "_", membership)
  
  CLL_data_psms_filt <- tryCatch(
    get_CLL_and_annotate(CLL_data_psms = cll_psms, graphs_comms = graphs_comms, ioi = gene_symbol),
    error = function(e) { data.frame() }  #empty data frame on error
  )
  
  modifications_list <- list()
  
  process_data <- function(data, mode) {
    modifications <- list()
    if (nrow(data) > 0) {
      if (mode == "gene_symbol") {
        data$membership <- 0
      }
      
      all_types <- unique(data$type)
      for (type in all_types) {
        type_specific_data <- data[data$type == type , ]
        if(mode != "gene_symbol"){
          
          type_specific_data <- type_specific_data[type_specific_data$membership == membership , ]
        }
        if("tmt_value" %in% names(type_specific_data) && !is.null(type_specific_data$tmt_value) && length(type_specific_data$tmt_value) > 0) {
          type_sum <- sum(type_specific_data$tmt_value, na.rm = TRUE)
          
          if (type_sum != 0) {
            modifications[[type]] <- type_sum
          }
        }
      }
    }
    return(modifications)
  }
  

  randomized_data <- randomize_data2(CLL_data_psms_filt)
  modifications_list[["random"]] <- process_data(randomized_data, "random")
  
  return(list(id = subset_id, modifications = modifications_list))
}


numCores <- detectCores() - 1
cl <- makeCluster(numCores)

#Run in parallel cluster
clusterExport(cl, list("process_subset2", "cll_psms","randomize_data2", "get_CLL_and_annotate", "graphs_comms", "all_from_CLL", "unique_ids"))
clusterEvalQ(cl, {
  library(dplyr)
  library(tidyr)
  library(igraph)
})


results <- parLapply(cl, unique_ids, function(id) {
  parts <- strsplit(id, "_", fixed = TRUE)[[1]]
  gene_symbol <- parts[1]
  membership <- parts[2]
  
  
  
  process_subset2(gene_symbol, membership, cll_psms, graphs_comms, all_from_CLL)
})

all_from_CLL_random2 <- all_from_CLL[0, ]
for (result in results) {
  id <- result$id
  modifications_list <- result$modifications
  
  all_from_CLL_random2 <- update_mode_df(all_from_CLL_random2, modifications_list[["random"]], id, all_from_CLL)
}

all_from_CLL_random2_normcenter <- normalize_and_center(all_from_CLL_random2)
stopCluster(cl)



saveRDS(all_from_CLL_random2_normcenter, "data/all_from_CLL_random2_normcenter.RDS")
