library(igraph)
library(tidyr)
library(dplyr)
library(vsn)
library(parallel)
#library(MultiAssayExperiment)
'%!in%' <- function(x,y)!('%in%' (x,y))
expand_rows_with_multiple_terms <- function(df, column_name) {
  
  df <- df %>%
    mutate(!!column_name := strsplit(as.character(.[[column_name]]), ",")) %>%
    unnest(!!column_name)
  
  return(df)
}
setwd("~/Documents/githubspot/ibrutinib-cll/")

graphs_comms <- readRDS("~/Documents/githubspot/ibrutinib-cll/data/graphs_comms_again.RDS")

hits_everything <- readRDS("~/Documents/githubspot/ibrutinib-cll/output/hits_everything.RDS")



cll_meta <- read.delim("~/Downloads/CLL_metadata 1.txt") 
#cll_metadata_shiny <- load("~/Downloads/multiomics_MAE.RData")
#instead, separately processed from MAE (to avoid dependency clashes)

cll_meta <- expand_rows_with_multiple_terms(cll_meta, "TMT")
all_from_CLL <- readRDS("data/all_from_CLL_default_normcenter.RDS")
all_from_CLL_normcenter <- readRDS("data/all_from_CLL_default_normcenter.RDS") #some code expects this name (merge from other code)

big_CLL_meta <- read.csv("CLL_clinical_sample_metadata.csv")
cll_translation_table <- read.delim("~/Downloads/CLL_patIDmeta_translation_table.txt")

big_CLL_meta <- left_join(big_CLL_meta, cll_translation_table, by = "patient_ID")
big_CLL_meta$TMT <- big_CLL_meta$TMT.channel
cll_meta <- left_join(cll_meta, big_CLL_meta, by = "TMT")

#treated with ibrutinib
cll_meta$treated <- ifelse(grepl("Ibrutinib", cll_meta$specification_lastpretreatment) | grepl("Ibrutinib", cll_meta$specification_current_treatment), "Ibrutinib", "Other")

#treatment_status
cll_meta <- cll_meta %>%
  mutate(
    treatment_status = case_when(
      grepl("Ibrutinib", specification_lastpretreatment) ~ "pretreated",
      grepl("Ibrutinib", specification_current_treatment) ~ "intreatment",
      TRUE ~ "other" 
    )
  )

# cll_meta <- cll_meta %>%
#   mutate(
#     treatment_status_number = case_when(
#       grepl(3, pretreatment) ~ "pretreated",
#       grepl(3, current_treatment) ~ "intreatment",
#       TRUE ~ "other" 
#     )
#   )
# View(cll_meta[,c("treatment_status_number", "patient_ID.x")] %>% distinct()) 
####################################################################################
####################################################################################









####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
cll_proteoforms <- readRDS("data/all_from_CLL_default_normcenter.RDS")
#subset + process
cll_proteoforms <- cll_proteoforms[!is.na(rowSums(as.matrix(select(cll_proteoforms, starts_with("Set"))), na.rm = TRUE)),] %>%
  select(starts_with(c("Set", "id", "gene.sym"))) %>%
  distinct()
all_from_CLL <- cll_proteoforms
#Plot normalized
sample_columns <- colnames(all_from_CLL)[grepl("Set", colnames(all_from_CLL))]
 all_from_CLL_long <- stack(all_from_CLL[sample_columns])
 all_from_CLL_long$values <- as.numeric(all_from_CLL_long$values)
library(ggplot2)
library(ggpubr)
library(ggsignif)

#name the columns in each group
rest_samples <- unique(cll_meta[cll_meta$treated == "Other",]$TMT)
ibrutinib_samples <- unique(cll_meta[cll_meta$treated == "Ibrutinib",]$TMT)

all_from_CLL_long$id <- all_from_CLL$id
all_from_CLL_long$gene.symbol <- all_from_CLL$gene.symbol
all_from_CLL_long$Treatment <- all_from_CLL_long$id %in% ibrutinib_samples


plot_samples <- function(id, ibrutinib_samples, rest_samples, all_from_CLL, mode = "Gene") {

  
  if (!id %in% all_from_CLL$id) {
    stop(paste("ID", id, "not found in the dataframe column names."))
  }
  specific_row <- all_from_CLL %>% filter(id == !!id)
  
  selected_columns <- c("id", ibrutinib_samples, rest_samples)
  data_long <- specific_row %>%
    select(all_of(selected_columns)) %>%
    pivot_longer(cols = -id, names_to = "Sample", values_to = "Value") %>%
    mutate(Condition = case_when(
      Sample %in% ibrutinib_samples ~ "Ibrutinib",
      Sample %in% rest_samples ~ "Rest"
    ))
  
  #Create the box plot
  p <- ggplot(data_long, aes(x = Condition, y = Value, fill = Condition)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = paste("Comparison for ID:", id),
         x = "Condition",
         y = "Value") +
    scale_fill_manual(values = c("Ibrutinib" = "skyblue", "Rest" = "lightgreen")) +
    stat_compare_means(method = "t.test", label = "p.signif", 
                         vjust = 0.5,      
                         tip.length = 0.02
    )
  
  

  print(p)
}
plot_treatment_comparisons <- function(all_from_CLL, treatment, id, 
                                       treated_factors_standard = "TRUE",
                                       highlight_pretreated = FALSE) {
  
  treatment <- treatment %>% distinct(ColumnName, .keep_all = TRUE)
  
  specific_row <- all_from_CLL %>% filter(id == !!id)
  
  data_long <- specific_row %>%
    select(id, all_of(treatment$ColumnName)) %>%
    pivot_longer(cols = -id, names_to = "ColumnName", values_to = "Value") %>%
    left_join(treatment, by = "ColumnName") %>%
    select(-id, -ColumnName)

  if(n_distinct(data_long$Condition) < 2) {
    stop("Error: not enough unique conditions for the provided ID to perform comparison.")
  }
  
  #create comparisons list
  my_comparisons <- combn(unique(treatment$Condition), 2)
  my_comparisons <- lapply(1:ncol(my_comparisons), function(i) my_comparisons[, i])
  
  if(treated_factors_standard == "TRUE"){
    data_long$Condition <- factor(data_long$Condition, levels = c("other", "intreatment", "pretreated"))
    my_comparisons <- list(my_comparisons[[1]], my_comparisons[[3]], my_comparisons[[2]])
    }
  #Create the box plot
  p <- ggplot(data_long, aes(x = Condition, y = Value, fill = Condition)) + 
    geom_boxplot() +
    theme_pubclean() +
    labs(title = paste(id), 
         x = "Treatment Condition", 
         y = "Value") +
    scale_fill_brewer(palette = "Paired") + 
    facet_grid(. ~ ., scales = "free") + 
    guides(fill = "none") +
    xlab(NULL) +
    ylab(NULL) +
    geom_signif(comparisons = my_comparisons, 
                  map_signif_level=TRUE,
                step_increase = 0.1,
                test= "wilcox.test")
  
  #Process highlight
  
  if(highlight_pretreated == TRUE)
    p <- p +geom_point(data = subset(data_long, highlight_pretreated != ""), color = "black", size = 10) +
    geom_point(data = subset(data_long, highlight_pretreated != ""), color = "white", size = 9) +
    geom_point(data = subset(data_long, highlight_pretreated != ""), color = "#A6CEE3", size = 5) 

    return(p)
}
generate_treatment_comparisons <- function(all_from_CLL, treatment, id, treated_factors_standard = TRUE) {
  
  treatment <- treatment %>% distinct(ColumnName, .keep_all = TRUE)
  
  specific_row <- all_from_CLL %>% filter(id == !!id)
  
  data_long <- specific_row %>%
    select(id, all_of(treatment$ColumnName)) %>%
    pivot_longer(cols = -id, names_to = "ColumnName", values_to = "Value") %>%
    left_join(treatment, by = "ColumnName") %>%
    select(-id, -ColumnName)
  
  if(n_distinct(data_long$Condition) < 2) {
    stop("Error: not enough unique conditions for the provided ID to perform comparison.")
  }
  
  my_comparisons <- combn(unique(treatment$Condition), 2)
  my_comparisons <- lapply(1:ncol(my_comparisons), function(i) my_comparisons[, i])
  
  #Create the box plot
  p <- pairwise.wilcox.test(data_long$Value, data_long$Condition,
                            p.adjust.method = "BH", alternative = "two.sided")
  
  return(p)
}



#Create wilcox test data frame
###
#Pretreated separate
treatment = data.frame("Condition" = cll_meta$treatment_status,"ColumnName" = cll_meta$TMT)
wilcoxtest_df <- distinct(data.frame("id" = all_from_CLL_normcenter$id,
                            "other.v.intreatment" = "N.D.",
                            "other.v.pretreated" = "N.D.",
                            "intreatment.v.pretreated" = "N.D."))

unique_ids <- unique(wilcoxtest_df$id)
for (i in unique_ids) {
  do_wilcoxtest <- tryCatch(
    generate_treatment_comparisons(all_from_CLL = distinct(all_from_CLL_normcenter), treatment, id = i, treated_factors_standard = TRUE),
    error = function(e) {
      cat(paste("Error for ID", i, ": ", conditionMessage(e), "\n"))
      return(NULL) 
    }
  )
  
  if (!is.null(do_wilcoxtest)) {
   
      wilcoxtest_df[wilcoxtest_df$id == i, "other.v.intreatment"] <- do_wilcoxtest$p.value[1, 1]
      wilcoxtest_df[wilcoxtest_df$id == i, "other.v.pretreated"] <- do_wilcoxtest$p.value[2, 1]
      wilcoxtest_df[wilcoxtest_df$id == i, "intreatment.v.pretreated"] <- do_wilcoxtest$p.value[2, 2]
    
  }
}

write.csv(wilcoxtest_df, "wilcoxtest_df.csv")

#Which was a hit by f-teest
f_test <- read.csv("f_test_data_significant_null.csv")
priority_hits <- unique(wilcoxtest_df$id[wilcoxtest_df$id %in% f_test$id&
                                    (wilcoxtest_df$other.v.intreatment<0.05|
                                    wilcoxtest_df$other.v.intreatment<0.05|
                                    wilcoxtest_df$other.v.intreatment<0.05)])
###
#Pretreated in treatment together
treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT)
wilcoxtest_df_all <- data.frame("id" = all_from_CLL_normcenter$id,
                            "other.v.all.ibrutinib" = "N.D.")

unique_ids <- unique(wilcoxtest_df_all$id)
for (i in unique_ids) {
  do_wilcoxtest <- tryCatch(
    generate_treatment_comparisons(all_from_CLL = distinct(all_from_CLL_normcenter), treatment, id = i, treated_factors_standard = FALSE),
    error = function(e) {
      cat(paste("Error for ID", i, ": ", conditionMessage(e), "\n"))
      return(NULL) 
    }
  )
  
  if (!is.null(do_wilcoxtest)) {
    
    wilcoxtest_df_all[wilcoxtest_df_all$id == i, "other.v.all.ibrutinib"] <- do_wilcoxtest$p.value[1, 1]
    
  }
}
wilcoxtest_df_all <- distinct(wilcoxtest_df_all)

write.csv(wilcoxtest_df_all, "wilcoxtest_df_all.csv")
#wilcoxtest_df_all$adjusted_p <- p.adjust(wilcoxtest_df_all$other.v.all.ibrutinib, method = "BH")
#Don't do this - the function has already BH adjusted

priority_hits_all <- unique(wilcoxtest_df_all$id[wilcoxtest_df_all$id %in% f_test$id &
                                           as.numeric(wilcoxtest_df_all$other.v.all.ibrutinib)<0.05])


###########
###############
#######################
#IMAGES FOR PLOTS

pdf("washc2c_1_CLL.pdf", width = 3.63, height =  3.56)
plot_treatment_comparisons(id = "WASHC2C_1", 
                           treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT, "highlight_pretreated" = 
                                                    ifelse(cll_meta$treatment_status != "pretreated", "", "pretreated")), 
                           all_from_CLL = all_from_CLL_normcenter, 
                           treated_factors_standard = F,
                           highlight_pretreated = T)
dev.off()

pdf("washc2c_2_CLL.pdf", width = 3.63, height =  3.56)
plot_treatment_comparisons(id = "WASHC2C_2", 
                           treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT, "highlight_pretreated" = 
                                                    ifelse(cll_meta$treatment_status != "pretreated", "", "pretreated")), 
                           all_from_CLL = all_from_CLL_normcenter, 
                           treated_factors_standard = F,
                           highlight_pretreated = T)
dev.off()
pdf("SCYL1_2_CLL.pdf", width = 3.63, height =  3.56)
plot_treatment_comparisons(id = "SCYL1_2", 
                           treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT, "highlight_pretreated" = 
                                                    ifelse(cll_meta$treatment_status != "pretreated", "", "pretreated")), 
                           all_from_CLL = all_from_CLL_normcenter, 
                           treated_factors_standard = F,
                           highlight_pretreated = T)

dev.off()
pdf("CARM1_2_CLL.pdf", width = 3.63, height =  3.56)
plot_treatment_comparisons(id = "CARM1_2", 
                           treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT, "highlight_pretreated" = 
                                                    ifelse(cll_meta$treatment_status != "pretreated", "", "pretreated")), 
                           all_from_CLL = all_from_CLL_normcenter, 
                           treated_factors_standard = F,
                           highlight_pretreated = T)

dev.off()

pdf("YWHAE_2_CLL.pdf", width = 3.63, height =  3.56)
plot_treatment_comparisons(id = "YWHAE_2", 
                           treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT, "highlight_pretreated" = 
                                                    ifelse(cll_meta$treatment_status != "pretreated", "", "pretreated")), 
                           all_from_CLL = all_from_CLL_normcenter, 
                           treated_factors_standard = F,
                           highlight_pretreated = T)

dev.off()

pdf("BRAF_2_CLL.pdf", width = 3.63, height =  3.56)
plot_treatment_comparisons(id = "BRAF_2", 
                           treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT, "highlight_pretreated" = 
                                                    ifelse(cll_meta$treatment_status != "pretreated", "", "pretreated")), 
                           all_from_CLL = all_from_CLL_normcenter, 
                           treated_factors_standard = F,
                           highlight_pretreated = T)

dev.off()
pdf("BRAF_1_CLL.pdf", width = 3.63, height =  3.56)
plot_treatment_comparisons(id = "BRAF_1", 
                           treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT, "highlight_pretreated" = 
                                                    ifelse(cll_meta$treatment_status != "pretreated", "", "pretreated")), 
                           all_from_CLL = all_from_CLL_normcenter, 
                           treated_factors_standard = F,
                           highlight_pretreated = T)

dev.off()
pdf("ARAF_1_CLL.pdf", width = 3.63, height =  3.56)
plot_treatment_comparisons(id = "ARAF_1", 
                           treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT, "highlight_pretreated" = 
                                                    ifelse(cll_meta$treatment_status != "pretreated", "", "pretreated")), 
                           all_from_CLL = all_from_CLL_normcenter, 
                           treated_factors_standard = F,
                           highlight_pretreated = T)

dev.off()
pdf("COG3_1_CLL.pdf", width = 4.63, height =  3.56)
plot_treatment_comparisons(id = "COG3_1", 
                           treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT, "highlight_pretreated" = 
                                                    ifelse(cll_meta$treatment_status != "pretreated", "", "pretreated")), 
                           all_from_CLL = all_from_CLL_normcenter, 
                           treated_factors_standard = F,
                           highlight_pretreated = T)
dev.off()

pdf("COG1_1_CLL.pdf", width = 4.63, height =  3.56)
plot_treatment_comparisons(id = "COG1_1", 
                           treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT, "highlight_pretreated" = 
                                                    ifelse(cll_meta$treatment_status != "pretreated", "", "pretreated")), 
                           all_from_CLL = all_from_CLL_normcenter, 
                           treated_factors_standard = F,
                           highlight_pretreated = T)
dev.off()
pdf("COG6_1_CLL.pdf", width = 4.63, height =  3.56)
plot_treatment_comparisons(id = "COG6_1", 
                           treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT, "highlight_pretreated" = 
                                                    ifelse(cll_meta$treatment_status != "pretreated", "", "pretreated")), 
                           all_from_CLL = all_from_CLL_normcenter, 
                           treated_factors_standard = F,
                           highlight_pretreated = T)
dev.off()


plot_treatment_comparisons(id = "COG3_1", treatment = data.frame("Condition" = cll_meta$treatment_status,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter)
plot_treatment_comparisons(id = "COG3_1", treatment = data.frame("Condition" = cll_meta$treatment_status,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter)
plot_treatment_comparisons(id = "BTK_1", treatment = data.frame("Condition" = cll_meta$treatment_status,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter)
plot_treatment_comparisons(id = "BTK_2", treatment = data.frame("Condition" = cll_meta$treatment_status,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter)
plot_treatment_comparisons(id = "WASHC2C_2", treatment = data.frame("Condition" = cll_meta$treatment_status,"ColumnName" = cll_meta$TMT), all_from_CLL = distinct(all_from_CLL_normcenter))
plot_treatment_comparisons(id = "WASHC2C_1", treatment = data.frame("Condition" = cll_meta$treatment_status,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter)
plot_treatment_comparisons(id = "BRAF_1", treatment = data.frame("Condition" = cll_meta$treatment_status,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter)
plot_treatment_comparisons(id = "BRAF_2", treatment = data.frame("Condition" = cll_meta$treatment_status,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter)

plot_treatment_comparisons(id = "BTK_2", treatment = data.frame("Condition" = cll_meta$treatment_status,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter)

#####
#Just two groups
plot_treatment_comparisons(id = "WASHC2C_2", treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter, treated_factors_standard = F)

plot_treatment_comparisons(id = "WASHC2C_1", treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter, treated_factors_standard = F)

plot_treatment_comparisons(id = "BTK_1", treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter, treated_factors_standard = F)
plot_treatment_comparisons(id = "BTK_2", treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter, treated_factors_standard = F)
plot_treatment_comparisons(id = "BRAF_2", treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter, treated_factors_standard = F)

plot_treatment_comparisons(id = "BRAF_1", treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter, treated_factors_standard = F)
plot_treatment_comparisons(id = "YES1_2", treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter, treated_factors_standard = F)

plot_treatment_comparisons(id = "YES1_1", treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter, treated_factors_standard = F)
plot_treatment_comparisons(id = "COG3_2", treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter, treated_factors_standard = F)

plot_treatment_comparisons(id = "COG3_1", treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter, treated_factors_standard = F)
plot_treatment_comparisons(id = "COG8_1", treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT), all_from_CLL = all_from_CLL_normcenter, treated_factors_standard = F)


#LOOK AT PG5 samples
pg5_samples <- unique(cll_meta[cll_meta$PG == "PG5",]$TMT)
non_pg5_samples <- unique(cll_meta[cll_meta$PG != "PG5",]$TMT)

#how many in each

treatment %>%
  group_by(Condition) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1L)

plot_treatment_comparisons(id = "BTK_1", 
                           treatment = data.frame("Condition" = paste0(cll_meta$treatment_status, ", PG5 = ", grepl("PG5", cll_meta$PG)),"ColumnName" = cll_meta$TMT), 
                           all_from_CLL = all_from_CLL_normcenter,
                           treated_factors_standard = FALSE)
plot_treatment_comparisons(id = "BTK_1", 
                           treatment = data.frame("Condition" = cll_meta$treatment_status,"ColumnName" = cll_meta$TMT)%>%
                             subset(., ColumnName %in% pg5_samples), 
                           all_from_CLL = all_from_CLL_normcenter[,colnames(all_from_CLL_normcenter) %in% c(pg5_samples, "id")])

plot_treatment_comparisons(id = "BTK_1", 
                           treatment = data.frame("Condition" = cll_meta$treatment_status,"ColumnName" = cll_meta$TMT)%>%
                             subset(., ColumnName %in% non_pg5_samples), 
                           all_from_CLL = all_from_CLL_normcenter[,colnames(all_from_CLL_normcenter) %in% c(non_pg5_samples, "id")])


####################################################################################
####################################################################################
#Correlate with coculture response to ibrutinib

ibrutinib_coculture_Herbst <- read.delim("data/ibrutinib_coculture_Herbst.txt")
ibrutinib_coculture_Herbst$value <- as.numeric(gsub(",", ".", ibrutinib_coculture_Herbst$value))
treatment = data.frame("Condition" = cll_meta$treated,"ColumnName" = cll_meta$TMT,
                       "CLL" = cll_meta$patient_ID.x)


generate_treatment_correlations <- function(all_from_CLL,type= "correlation", treatment, id, treated_factors_standard = TRUE, ibrutinib_coculture_Herbst) {
  
  treatment <- treatment %>% distinct(ColumnName, .keep_all = TRUE)
  
  specific_row <- all_from_CLL %>% filter(id == !!id)
  
  data_long <- specific_row %>%
    select(id, all_of(treatment$ColumnName)) %>%
    pivot_longer(cols = -id, names_to = "ColumnName", values_to = "Value") %>%
    left_join(treatment, by = "ColumnName") %>%
    select(-id, -ColumnName) %>% unique() %>% as.data.frame()
  data_long <- data_long[!duplicated(data_long$CLL), ]
  rownames(data_long) <- data_long$CLL
  
  #order in same order as drug screen
  data_long <- data_long[ibrutinib_coculture_Herbst$primary,]
  
  if(n_distinct(data_long$Condition) < 2) {
    stop("Error: not enough unique conditions")
  }

  #Correlation
  p <- cor(x = data_long$Value, y = ibrutinib_coculture_Herbst$value,
                            method = "spearman")
  
if(type == "full"){
  p <- cor.test(x = data_long$Value, y = ibrutinib_coculture_Herbst$value,
           method = "spearman")
  
}
  return(p)
}

generate_treatment_correlations(all_from_CLL = all_from_CLL_normcenter, treatment, id = "WASHC2C_1", ibrutinib_coculture_Herbst = ibrutinib_coculture_Herbst)

spearmancorr_df <- data.frame("id" = all_from_CLL_normcenter$id)

unique_ids <- unique(spearmancorr_df$id)
for (i in unique_ids) {
  do_spearmancorr <- tryCatch(
    generate_treatment_correlations(all_from_CLL = all_from_CLL_normcenter, 
                                    treatment, id = i, 
                                    treated_factors_standard = TRUE,
                                    ibrutinib_coculture_Herbst = ibrutinib_coculture_Herbst),
    error = function(e) {
      cat(paste("Error for ID", i, ": ", conditionMessage(e), "\n"))
      return(NULL) 
    }
  )
  
  if (!is.null(do_spearmancorr)) {
    
    spearmancorr_df[spearmancorr_df$id == i, "all"] <- generate_treatment_correlations(all_from_CLL = all_from_CLL_normcenter, 
                                                                                       treatment, id = i, 
                                                                                       treated_factors_standard = TRUE,
                                                                                       ibrutinib_coculture_Herbst = ibrutinib_coculture_Herbst)
    spearmancorr_df[spearmancorr_df$id == i, "pg5"] <- generate_treatment_correlations(all_from_CLL = 
                                                                                         all_from_CLL_normcenter[,colnames(all_from_CLL_normcenter) %in% c(pg5_samples, "id")], 
                                                                                             treatment%>%subset(., ColumnName %in% pg5_samples), id = i, 
                                                                                             treated_factors_standard = TRUE,
                                                                                             ibrutinib_coculture_Herbst = ibrutinib_coculture_Herbst[ibrutinib_coculture_Herbst$PG == 5,])
    spearmancorr_df[spearmancorr_df$id == i, "nonpg5"] <- generate_treatment_correlations(all_from_CLL = 
                                                                                            all_from_CLL_normcenter[,colnames(all_from_CLL_normcenter) %in% c(non_pg5_samples, "id")], 
                                                                                              treatment%>%
                                                                                            subset(., ColumnName %in% non_pg5_samples), id = i, 
                                                                                              treated_factors_standard = TRUE,
                                                                                              ibrutinib_coculture_Herbst = ibrutinib_coculture_Herbst[ibrutinib_coculture_Herbst$PG != 5,])
    
  }
}

spearmancorr_df <- distinct(spearmancorr_df)
saveRDS(spearmancorr_df, file = "data/spearmancorr_df.RDS")

compare_correlations <- function(r1, r2, n1, n2) {
  #Fisher Z transformation
  Z1 <- 0.5 * log((1 + r1) / (1 - r1))
  Z2 <- 0.5 * log((1 + r2) / (1 - r2))
  
  #Standard error of the difference
  SE_difference <- sqrt(1/(n1-3) + 1/(n2-3))
  
  #Z-score for the difference
  Z_score <- (Z1 - Z2) / SE_difference
  
  #P-value from the Z-score
  p_value <- 2 * (1 - pnorm(abs(Z_score)))
  
  return(list(Z_score = Z_score, P_value = p_value))
}

spearmancorr_df$comparsion_corr_Z <- NA
spearmancorr_df$comparsion_corr_P <- NA

for(i in spearmancorr_df$id){
  
  spearmancorr_df[spearmancorr_df$id == i,]$comparsion_corr_Z <- compare_correlations(spearmancorr_df[spearmancorr_df$id == i,]$nonpg5,
                                                                                    spearmancorr_df[spearmancorr_df$id == i,]$pg5,
                                                                                    length(non_pg5_samples),length(pg5_samples))$Z_score

  
  spearmancorr_df[spearmancorr_df$id == i,]$comparsion_corr_P <- compare_correlations(spearmancorr_df[spearmancorr_df$id == i,]$nonpg5,
                                                                                      spearmancorr_df[spearmancorr_df$id == i,]$pg5,
                                                                                      length(non_pg5_samples),length(pg5_samples))$P_value 
}


spearmancorr_df$neg_log_p <- -log10(spearmancorr_df$comparsion_corr_P)
significance_threshold <- -log10(0.05) 
z_score_threshold <- 1.96

#Create volcano plot
ggplot(spearmancorr_df, aes(x = comparsion_corr_Z, y = neg_log_p)) +
  geom_point(aes(color = neg_log_p), size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Z-score",
       y = "-log10(P-value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = significance_threshold, linetype = "dashed", color = "darkgrey") +  #Add a reference line for P-value = 0.05
  geom_vline(xintercept = c(-z_score_threshold, z_score_threshold), linetype = "dashed", color = "darkgrey") +  #Add reference lines for significance thresholds
  ggrepel::geom_label_repel(data = subset(spearmancorr_df, neg_log_p > significance_threshold & abs(comparsion_corr_Z) > z_score_threshold),
            aes(label = id), max.overlaps = 10) +
  geom_point(data = subset(spearmancorr_df, id %in% hits_everything[hits_everything$pAdj < 0.05,]$id), color = "black", size = 10) +
  geom_point(data = subset(spearmancorr_df, id %in% hits_everything[hits_everything$pAdj < 0.05,]$id), color = "white", size = 9) +
  geom_point(data = subset(spearmancorr_df, id %in% hits_everything[hits_everything$pAdj < 0.05,]$id), color = "#A6CEE3", size = 5) 


ggplot(subset(spearmancorr_df, id %in% hits_everything[hits_everything$pAdj < 0.05,]$id), aes(x = comparsion_corr_Z, y = neg_log_p)) +
  geom_point(aes(color = neg_log_p), size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Z-score",
       y = "-log10(P-value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = significance_threshold, linetype = "dashed", color = "darkgrey") +  #Add a reference line for P-value = 0.05
  geom_vline(xintercept = c(-z_score_threshold, z_score_threshold), linetype = "dashed", color = "darkgrey") +  #Add reference lines for significance thresholds
  ggrepel::geom_label_repel(data = subset(subset(spearmancorr_df, id %in% hits_everything[hits_everything$pAdj < 0.05,]$id), neg_log_p > significance_threshold & abs(comparsion_corr_Z) > z_score_threshold),
                            aes(label = id), max.overlaps = 10) +
  geom_point(data = subset(spearmancorr_df, id %in% hits_everything[hits_everything$pAdj < 0.05,]$id), color = "black", size = 10) +
  geom_point(data = subset(spearmancorr_df, id %in% hits_everything[hits_everything$pAdj < 0.05,]$id), color = "white", size = 9) +
  geom_point(data = subset(spearmancorr_df, id %in% hits_everything[hits_everything$pAdj < 0.05,]$id), color = "#A6CEE3", size = 5) 


#Plot correlations
rownames(cll_proteoforms) <- cll_proteoforms$id
data_long <- cll_proteoforms %>%
  select(id, all_of(treatment$ColumnName)) %>%
  pivot_longer(cols = -id, names_to = "ColumnName", values_to = "Value") %>%
  left_join(treatment, by = "ColumnName") %>%
  unique() %>% as.data.frame()

merged_data <- merge(data_long, ibrutinib_coculture_Herbst, by.x = "CLL", by.y = "primary")

i = "BZW2_2"
merged_data_i <- merged_data[merged_data$id == i, ]
merged_data_i$Value <- as.numeric(merged_data_i$Value)
merged_data_i$value <- as.numeric(merged_data_i$value)
merged_data_i$PG <- as.factor(merged_data_i$PG)

#Create scatterplot
ggplot(merged_data_i, aes(x = Value, y = value)) +
  geom_point(aes(color = PG5), size = 3) +
  scale_color_brewer(palette = "Set1")+
  geom_smooth(method = "lm", color = "grey", se = TRUE) + 
  geom_smooth(method = "lm", color = "red", se = TRUE, data = subset(merged_data_i, ColumnName %in% pg5_samples)) + 
  labs(title = paste("PG5 Spearman Correlation: ", round(spearmancorr_df[spearmancorr_df$id == i ,]$pg5, 2)),
       x = "Proteoform abundance",
       y = "Ibrutinib coculture") +
  theme_minimal() +
  theme(legend.position = "none")# + 
  #geom_label(aes(label=PG5), vjust = -1, hjust=1.1, size=3)

i = "BZW2_2"
merged_data_i <- merged_data[merged_data$id == i, ]
merged_data_i$Value <- as.numeric(merged_data_i$Value)
merged_data_i$value <- as.numeric(merged_data_i$value)
merged_data_i$PG <- as.factor(merged_data_i$PG)

#Create scatterplot
pdf(width = 2.43, height = 3.09, "ASB_rest_BZW2_2.pdf")
ggplot(merged_data_i, aes(x = Value, y = value)) +
  geom_point(aes(color = PG5), size = 3) +
  scale_color_brewer(palette = "Set1")+
  geom_smooth(method = "lm", color = "grey", se = TRUE) + 
  geom_smooth(method = "lm", color = "red", se = TRUE, data = subset(merged_data_i, ColumnName %in% pg5_samples)) + 
  labs(#title = paste("High-Splice Spearman Correlation: ", round(spearmancorr_df[spearmancorr_df$id == i ,]$pg5, 2)),
       x = "Proteoform abundance",
       y = "Ibrutinib coculture") +
  theme_minimal() +
  theme(legend.position = "none")# + 
#geom_label(aes(label=PG5), vjust = -1, hjust=1.1, size=3)
dev.off()
#separate
pdf(width = 2.43, height = 3.09, "ASB_BZW2_2.pdf")
ggplot(merged_data_i[merged_data_i$PG5 == "ASB-CLL",], aes(x = Value, y = value)) +
  geom_point(aes(color = PG5), size = 3) +
  scale_color_brewer(palette = "Set1")+
  #geom_smooth(method = "lm", color = "grey", se = TRUE) + 
  geom_smooth(method = "lm", color = "red", se = TRUE, data = subset(merged_data_i, ColumnName %in% pg5_samples)) + 
  labs(#title = paste("High-Splice Spearman Correlation: ", round(spearmancorr_df[spearmancorr_df$id == i ,]$pg5, 2)),
    x = "Proteoform abundance",
    y = "Ibrutinib coculture") +
  theme_minimal() +
  theme(legend.position = "none")# + 
#geom_label(aes(label=PG5), vjust = -1, hjust=1.1, size=3)
dev.off()

pdf(width = 2.43, height = 3.09, "rest_BZW2_2.pdf")
ggplot(merged_data_i[merged_data_i$PG5 != "ASB-CLL",], aes(x = Value, y = value)) +
  geom_point(aes(color = PG5), size = 3, color = "#161FFF") +
  geom_smooth(method = "lm", color = "grey", se = TRUE) + 
  #geom_smooth(method = "lm", color = "red", se = TRUE, data = subset(merged_data_i, ColumnName %in% pg5_samples)) + 
  labs(#title = paste("High-Splice Spearman Correlation: ", round(spearmancorr_df[spearmancorr_df$id == i ,]$pg5, 2)),
    x = "Proteoform abundance",
    y = "Ibrutinib coculture") +
  theme_minimal() +
  theme(legend.position = "none")# + 
#geom_label(aes(label=PG5), vjust = -1, hjust=1.1, size=3)
dev.off()

i = "BZW2_1"
merged_data_i <- merged_data[merged_data$id == i, ]
merged_data_i$Value <- as.numeric(merged_data_i$Value)
merged_data_i$value <- as.numeric(merged_data_i$value)
merged_data_i$PG <- as.factor(merged_data_i$PG)

#Create scatterplot
pdf(width = 2.43, height = 3.09, "ASB_rest_BZW2_1.pdf")
ggplot(merged_data_i, aes(x = Value, y = value)) +
  geom_point(aes(color = PG5), size = 3) +
  scale_color_brewer(palette = "Set1")+
  geom_smooth(method = "lm", color = "grey", se = TRUE) + 
  geom_smooth(method = "lm", color = "red", se = TRUE, data = subset(merged_data_i, ColumnName %in% pg5_samples)) + 
  labs(#title = paste("High-Splice Spearman Correlation: ", round(spearmancorr_df[spearmancorr_df$id == i ,]$pg5, 2)),
    x = "Proteoform abundance",
    y = "Ibrutinib coculture") +
  theme_minimal() +
  theme(legend.position = "none")# + 
#geom_label(aes(label=PG5), vjust = -1, hjust=1.1, size=3)
dev.off()
#separate
pdf(width = 2.43, height = 3.09, "ASB_BZW2_1.pdf")
ggplot(merged_data_i[merged_data_i$PG5 == "ASB-CLL",], aes(x = Value, y = value)) +
  geom_point(aes(color = PG5), size = 3) +
  scale_color_brewer(palette = "Set1")+
  #geom_smooth(method = "lm", color = "grey", se = TRUE) + 
  geom_smooth(method = "lm", color = "red", se = TRUE, data = subset(merged_data_i, ColumnName %in% pg5_samples)) + 
  labs(#title = paste("High-Splice Spearman Correlation: ", round(spearmancorr_df[spearmancorr_df$id == i ,]$pg5, 2)),
    x = "Proteoform abundance",
    y = "Ibrutinib coculture") +
  theme_minimal() +
  theme(legend.position = "none")# + 
#geom_label(aes(label=PG5), vjust = -1, hjust=1.1, size=3)
dev.off()

pdf(width = 2.43, height = 3.09, "rest_BZW2_1.pdf")
ggplot(merged_data_i[merged_data_i$PG5 != "ASB-CLL",], aes(x = Value, y = value)) +
  geom_point(aes(color = PG5), size = 3, color = "#161FFF") +
  geom_smooth(method = "lm", color = "grey", se = TRUE) + 
  #geom_smooth(method = "lm", color = "red", se = TRUE, data = subset(merged_data_i, ColumnName %in% pg5_samples)) + 
  labs(#title = paste("High-Splice Spearman Correlation: ", round(spearmancorr_df[spearmancorr_df$id == i ,]$pg5, 2)),
    x = "Proteoform abundance",
    y = "Ibrutinib coculture") +
  theme_minimal() +
  theme(legend.position = "none")# + 
#geom_label(aes(label=PG5), vjust = -1, hjust=1.1, size=3)
dev.off()
merged_data_i <- merged_data_i[order(merged_data_i$PG5, decreasing = FALSE),]
pdf(width = 3.43, height = 3.09, "ASB_CLL_ibrutinib_result.pdf")
ggplot(merged_data_i, aes(x = PG5, y = value)) +
  geom_boxplot(aes(fill = PG5)) +
  ggthemes::scale_fill_canva()+
  theme_minimal() +
  geom_hline(aes(yintercept =1, color = "grey")) +
  theme(legend.position = "none") + 
  xlab("")
dev.off()
spearmancorr_df_f <- spearmancorr_df[spearmancorr_df$id %in% f_test$id,]

i = "EIF2B5_1"
merged_data_i <- merged_data[merged_data$id == i, ]
merged_data_i$Value <- as.numeric(merged_data_i$Value)
merged_data_i$value <- as.numeric(merged_data_i$value)
merged_data_i$PG <- as.factor(merged_data_i$PG)

#Create scatterplot
pdf(width = 2.43, height = 3.09, "ASB_rest_EIF2B5_1.pdf")
ggplot(merged_data_i, aes(x = Value, y = value)) +
  geom_point(aes(color = PG5), size = 3) +
  scale_color_brewer(palette = "Set1")+
  geom_smooth(method = "lm", color = "grey", se = TRUE) + 
  geom_smooth(method = "lm", color = "red", se = TRUE, data = subset(merged_data_i, ColumnName %in% pg5_samples)) + 
  labs(#title = paste("High-Splice Spearman Correlation: ", round(spearmancorr_df[spearmancorr_df$id == i ,]$pg5, 2)),
    x = "Proteoform abundance",
    y = "Ibrutinib coculture") +
  theme_minimal() +
  theme(legend.position = "none")# + 
#geom_label(aes(label=PG5), vjust = -1, hjust=1.1, size=3)
dev.off()

i = "EIF2B4_1"
merged_data_i <- merged_data[merged_data$id == i, ]
merged_data_i$Value <- as.numeric(merged_data_i$Value)
merged_data_i$value <- as.numeric(merged_data_i$value)
merged_data_i$PG <- as.factor(merged_data_i$PG)

#Create scatterplot
pdf(width = 2.43, height = 3.09, "ASB_rest_EIF2B4_1.pdf")
ggplot(merged_data_i, aes(x = Value, y = value)) +
  geom_point(aes(color = PG5), size = 3) +
  scale_color_brewer(palette = "Set1")+
  geom_smooth(method = "lm", color = "grey", se = TRUE) + 
  geom_smooth(method = "lm", color = "red", se = TRUE, data = subset(merged_data_i, ColumnName %in% pg5_samples)) + 
  labs(#title = paste("High-Splice Spearman Correlation: ", round(spearmancorr_df[spearmancorr_df$id == i ,]$pg5, 2)),
    x = "Proteoform abundance",
    y = "Ibrutinib coculture") +
  theme_minimal() +
  theme(legend.position = "none")# + 
#geom_label(aes(label=PG5), vjust = -1, hjust=1.1, size=3)
dev.off()

i = "EIF5B_1"
merged_data_i <- merged_data[merged_data$id == i, ]
merged_data_i$Value <- as.numeric(merged_data_i$Value)
merged_data_i$value <- as.numeric(merged_data_i$value)
merged_data_i$PG <- as.factor(merged_data_i$PG)

#Create scatterplot
pdf(width = 2.43, height = 3.09, "ASB_rest_EIF5B1_1.pdf")
ggplot(merged_data_i, aes(x = Value, y = value)) +
  geom_point(aes(color = PG5), size = 3) +
  scale_color_brewer(palette = "Set1")+
  geom_smooth(method = "lm", color = "grey", se = TRUE) + 
  geom_smooth(method = "lm", color = "red", se = TRUE, data = subset(merged_data_i, ColumnName %in% pg5_samples)) + 
  labs(#title = paste("High-Splice Spearman Correlation: ", round(spearmancorr_df[spearmancorr_df$id == i ,]$pg5, 2)),
    x = "Proteoform abundance",
    y = "Ibrutinib coculture") +
  theme_minimal() +
  theme(legend.position = "none")# + 
#geom_label(aes(label=PG5), vjust = -1, hjust=1.1, size=3)
dev.off()

#Get the correlation p values
i = "BWZ2_2"

generate_treatment_correlations(all_from_CLL = all_from_CLL, 
                                                                                   treatment, id = i, 
                                                                                   treated_factors_standard = TRUE,type = "full",
                                                                                   ibrutinib_coculture_Herbst = ibrutinib_coculture_Herbst)
 generate_treatment_correlations(all_from_CLL = 
                                                                                     all_from_CLL_normcenter[,colnames(all_from_CLL_normcenter) %in% c(pg5_samples, "id")], 
                                                                                   treatment%>%subset(., ColumnName %in% pg5_samples), id = i, 
                                                                                   treated_factors_standard = TRUE, type = "full",
                                                                                   ibrutinib_coculture_Herbst = ibrutinib_coculture_Herbst[ibrutinib_coculture_Herbst$PG == 5,])
generate_treatment_correlations(all_from_CLL = 
                                                                                        all_from_CLL_normcenter[,colnames(all_from_CLL_normcenter) %in% c(non_pg5_samples, "id")], 
                                                                                      treatment%>%
                                                                                        subset(., ColumnName %in% non_pg5_samples), id = i, 
                                                                                      treated_factors_standard = TRUE,type = "full",
                                                                                      ibrutinib_coculture_Herbst = ibrutinib_coculture_Herbst[ibrutinib_coculture_Herbst$PG != 5,])



