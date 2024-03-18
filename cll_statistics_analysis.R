#library(igraph)
library(tidyr)
library(dplyr)
library(broom)
library(purrr)
library(ggplot2)
library(ggpubr)
#library(vsn)
#library(parallel)
library(readxl)
'%!in%' <- function(x,y)!('%in%' (x,y))
expand_rows_with_multiple_terms <- function(df, column_name) {
  
  df <- df %>%
    mutate(!!column_name := strsplit(as.character(.[[column_name]]), ",")) %>%
    unnest(!!column_name)
  
  return(df)
}
setwd("~/Documents/githubspot/ibrutinib-cll/")

#read in CLL proteforms
cll_proteoforms <- readRDS("data/all_from_CLL_default_normcenter.RDS")
    #subset + process
    cll_proteoforms <- cll_proteoforms[!is.na(rowSums(as.matrix(select(cll_proteoforms, starts_with("Set"))), na.rm = TRUE)),] %>%
      select(starts_with(c("Set", "id", "gene.sym"))) %>%
      distinct()

#read in CLL proteins
    #Using assembled from PSM, identically processed table
    cll_proteins <- readRDS("data/all_from_CLL_gene_symbol_normcenter.RDS")
    cll_proteins <- cll_proteins[!is.na(rowSums(as.matrix(select(cll_proteins, starts_with("Set"))), na.rm = TRUE)),] %>%
      select(starts_with(c("Set", "id", "gene.sym"))) %>%
      distinct()
    
    cll_null <- readRDS("data/all_from_CLL_random2_normcenter.RDS")
    
    cll_null <- cll_null[!is.na(rowSums(as.matrix(select(cll_proteins, starts_with("Set"))), na.rm = TRUE)),] %>%
      select(starts_with(c("Set", "id", "gene.sym"))) %>%
      distinct()
    #Commented code - using published proteins table
# cll_proteins <- readxl::read_xlsx("~/Downloads/41467_2022_33385_MOESM4_ESM.xlsx", sheet = 2) 
#     #subset + process
#     cll_proteins <- as.data.frame(cll_proteins[cll_proteins$`HGNC gene symbol` %in% cll_proteoforms$gene.symbol,])
#     rownames(cll_proteins) <- cll_proteins$`HGNC gene symbol`
#     cll_proteins <- cll_proteins %>% select(starts_with(c("CLL")))
#     
#     ####
#     #merge
#     cll_proteoforms <- cll_proteoforms[cll_proteoforms$gene.symbol %in% rownames(cll_proteins) & 
#                                          cll_proteoforms$id != "N.D.",]
#     cll_proteins_rownames <- rownames(cll_proteins)
#     rownames(cll_proteoforms) <- cll_proteoforms$id
#     genes_per_proteoform <- cll_proteoforms$gene.symbol
#     cll_proteoforms <- cll_proteoforms %>% 
#       select(starts_with(c("Set"))) 
#     cll_proteins_anno <- readxl::read_xlsx("~/Downloads/41467_2022_33385_MOESM4_ESM.xlsx", sheet = 3) 
#     colnames(cll_proteins) <- cll_proteins_anno$`TMT-channel`
        ####
        #merge
        cll_proteoforms <- cll_proteoforms[cll_proteoforms$gene.symbol %in% cll_proteins$gene.symbol &
                                             cll_proteoforms$id != "N.D.",]
        cll_proteins_rownames <- cll_proteins$gene.symbol
        rownames(cll_proteoforms) <- cll_proteoforms$id
        genes_per_proteoform <- cll_proteoforms$gene.symbol
        cll_proteoforms <- cll_proteoforms %>%
          select(starts_with(c("Set")))


#read in CLL meta
cll_meta <- read.delim("~/Downloads/CLL_metadata 1.txt") 
cll_meta <- expand_rows_with_multiple_terms(cll_meta, "TMT")


#read in drug screen
ibrutinib_coculture_Herbst <- read.delim("data/ibrutinib_coculture_Herbst.txt")
ibrutinib_coculture_Herbst$value <- as.numeric(gsub(",", ".", ibrutinib_coculture_Herbst$value))

#make analysis objects
treatment = data.frame("Condition" = cll_meta$treatment_status,"ColumnName" = cll_meta$TMT,
                       "CLL" = cll_meta$patient_ID)

#############################################
#ftest, protein versus proteoform

#variance for each proteoform
cll_proteoforms <- cll_proteoforms %>%
  rowwise() %>%
  mutate(variance_proteoform = var(c_across(starts_with("Set")), na.rm = TRUE)) %>%
  ungroup()

#add gene symbols
cll_proteoforms$gene_symbol <- genes_per_proteoform

#find maximum variance for each gene symbol at the proteoform level
max_variance_proteoforms <- cll_proteoforms %>%
  group_by(gene_symbol) %>%
  summarize(max_variance = if(all(is.na(variance_proteoform))) NA_real_ else max(variance_proteoform, na.rm = TRUE),
            .groups = 'drop')

#calculate variance for each gene symbol
cll_proteins <- cll_proteins %>%
  rowwise() %>%
  mutate(variance_gene_symbol = var(c_across(starts_with("Set")), na.rm = TRUE)) %>%
  ungroup()
cll_proteins$gene_symbol <- cll_proteins_rownames


#F-tests
max_variance_proteoforms <- max_variance_proteoforms %>%
  filter(!is.na(max_variance))
cll_proteins <- cll_proteins %>%
  filter(!is.na(variance_gene_symbol))

f_test_data <- left_join(max_variance_proteoforms, cll_proteins, by = "gene_symbol")

f_test_data <- f_test_data %>%
  mutate(f_statistic = max_variance / variance_gene_symbol,
         p_value = purrr::map2_dbl( max_variance, variance_gene_symbol, ~pf(.x / .y, df1 = ncol(cll_proteoforms) - 1, df2 = ncol(cll_proteins) - 1, lower.tail = FALSE)))

#reverse direction: compare gene versus proteoform
f_test_data <- f_test_data %>%
  mutate(f_statistic_reverse = variance_gene_symbol / max_variance,
         p_value_reverse = map2_dbl(variance_gene_symbol, max_variance, 
                                    ~pf(.x / .y, df1 = ncol(cll_proteins) - 1, 
                                        df2 = ncol(cll_proteoforms) - 1, lower.tail = FALSE)))
wilcox_test_result <- wilcox.test(f_test_data$f_statistic, 
                                  f_test_data$f_statistic_reverse, 
                                  paired = FALSE)


#Compare variances
f_test_data_long <- tidyr::pivot_longer(f_test_data,
                                        cols = c("max_variance", "variance_gene_symbol"), 
                                        names_to = "variance", 
                                        values_to = "variance_value")
f_test_data_long$variance <- gsub("variance_gene_symbol", "protein", f_test_data_long$variance )
f_test_data_long$variance <- gsub("max_variance", "proteoform", f_test_data_long$variance )
boxplot <- ggplot(f_test_data_long, aes(x = variance, y = variance_value, fill = variance)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#088DA5", "#E69F00")) +
  labs(title = "Variance Distributions", x = "Variance Type", y = "Variance") +
  theme_minimal() +
  #stat_compare_means(method = "wilcox.test", label = "p.signif") +
#  scale_y_log10()+
  theme(legend.position = "none") +
  guides(size = "none", alpha = "none", fill = "none")+
  theme_pubclean() +
  geom_signif(comparisons = list(c("protein", "proteoform")), 
              map_signif_level=FALSE,
              step_increase = 0.1,
              test= "wilcox.test")

#Print the plot
boxplot


#adjust p val
f_test_data$adjusted_p_values <- p.adjust(f_test_data$p_value, method = "BH")
f_test_data_protein <- f_test_data %>% select(., starts_with(c("gene_sym", "f_stat", "p_val", "adjusted_")))

#########################################
#ftest, proteoform random versus proteoform
#READ IN AGAIN
#read in CLL proteforms
cll_proteoforms <- readRDS("data/all_from_CLL_default_normcenter.RDS")
#subset + process
cll_proteoforms <- cll_proteoforms[!is.na(rowSums(as.matrix(select(cll_proteoforms, starts_with("Set"))), na.rm = TRUE)),] %>%
  select(starts_with(c("Set", "id", "gene.sym"))) %>%
  distinct()
rownames(cll_proteoforms) <- cll_proteoforms$id

cll_null <- readRDS("data/all_from_CLL_random2_normcenter.RDS")

cll_null <- cll_null[!is.na(rowSums(as.matrix(select(cll_null, starts_with("Set"))), na.rm = TRUE)),] %>%
  select(starts_with(c("Set", "id", "gene.sym"))) %>%
  distinct()
rownames(cll_null) <- cll_null$id
#variance for each proteoform
cll_proteoforms <- cll_proteoforms %>%
  rowwise() %>%
  mutate(variance_proteoform = var(c_across(starts_with("Set")), na.rm = TRUE)) %>%
  ungroup()

# #find maximum variance for each gene symbol at the proteoform level
# max_variance_proteoforms <- cll_proteoforms %>%
#   group_by(gene_symbol) %>%
#   summarize(max_variance = if(all(is.na(variance_proteoform))) NA_real_ else max(variance_proteoform, na.rm = TRUE),
#             .groups = 'drop')

#calculate variance for each gene symbol
cll_null <- cll_null %>%
  rowwise() %>%
  mutate(variance_random = var(c_across(starts_with("Set")), na.rm = TRUE)) %>%
  ungroup()


#F-tests
cll_null <- cll_null %>%
  filter(!is.na(variance_random))
cll_proteoforms <- cll_proteoforms %>%
  filter(!is.na(variance_proteoform))

f_test_data <- left_join(cll_proteoforms, cll_null, by = "id")

f_test_data <- f_test_data %>%
  mutate(f_statistic = variance_proteoform / variance_random,
         p_value = purrr::map2_dbl( variance_proteoform, variance_random, ~pf(.x / .y, df1 = ncol(cll_proteoforms) - 1, df2 = ncol(cll_null) - 1, lower.tail = FALSE)))

#Compare variances
f_test_data_long <- tidyr::pivot_longer(f_test_data,
                                        cols = c("variance_proteoform", "variance_random"), 
                                        names_to = "variance", 
                                        values_to = "variance_value")
f_test_data_long$variance <- gsub("variance_random", "random", f_test_data_long$variance )
f_test_data_long$variance <- gsub("variance_proteoform", "proteoform", f_test_data_long$variance )
boxplot <- ggplot(f_test_data_long, aes(x = variance, y = variance_value, fill = variance)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#088DA5", "#E69F00")) +
  labs(title = "Variance Distributions", x = "Variance Type", y = "Variance") +
  theme_minimal() +
  #stat_compare_means(method = "wilcox.test", label = "p.signif") +
  #  scale_y_log10()+
  theme(legend.position = "none") +
  guides(size = "none", alpha = "none", fill = "none")+
  theme_pubclean() +
  geom_signif(comparisons = list(c("random", "proteoform")), 
              map_signif_level=FALSE,
              step_increase = 0.1,
              test= "wilcox.test")

#Print the plot
boxplot

f_test_data <- f_test_data %>%
  mutate(f_statistic = variance_proteoform / variance_random,
         p_value = purrr::map2_dbl(variance_proteoform, variance_random, 
                                   ~pf(.x / .y, df1 = ncol(cll_proteoforms) - 1, 
                                       df2 = ncol(cll_null) - 1, lower.tail = FALSE))) %>%
  mutate(adjusted_p_values = p.adjust(p_value, method = "BH"),
         Significant = ifelse(adjusted_p_values < 0.05, "Yes", "No")) # Flag significant comparisons

#Prepare data for plotting
plot_data <- f_test_data %>%
  select(id, variance_proteoform, variance_random, Significant) %>%
  pivot_longer(cols = c(variance_proteoform, variance_random), 
               names_to = "Group", 
               values_to = "Variance")

plot_data$Group <- factor(plot_data$Group, levels = c("variance_proteoform", "variance_random"),
                          labels = c("Proteoform", "Random"))

plot <- ggplot(plot_data, aes(x = Group, y = Variance, color = Group)) +
  geom_boxplot() +
  geom_jitter(data = filter(plot_data, Significant == "Yes"), aes(color = Group, shape = Significant), 
              width = 0.2, size = 3, alpha = 1) +
  scale_shape_manual(values = c("Yes" = 10)) + 
  scale_color_manual(values = c("Proteoform" = "#088DA5", "Random" = "#E69F00")) +
  labs(title = "Variance Comparison: Proteoform vs. Random",
       x = "", y = "Variance") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  theme_pubclean() +
  geom_signif(comparisons = list(c("random", "proteoform")), 
              map_signif_level=FALSE,
              step_increase = 0.1,
              test= "wilcox.test")

#Print the plot
print(plot)

#plot
genesymbol_hits <- readRDS("data/re_search_all_from_CLL.RDS")
genesymbol_hits <-genesymbol_hits[genesymbol_hits$pAdj < 0.05,]
genesymbol_hits <- unique(genesymbol_hits$id[genesymbol_hits$id %in% f_test_data$id])

fstat_hits <- readRDS("shiny_stable6/nparc_df_rchacv.RDS")
fstat_hits <- fstat_hits %>%
  group_by(id) %>%
  slice(which.max(fStat)) %>%
  ungroup()

# fstat_hits$gene_symbol <- 0
# for (i in 1:nrow(fstat_hits)) {
#   parts <- strsplit(fstat_hits$id[i], "_", fixed = TRUE)[[1]]
#   if (!is.na(parts[1])) {
#     fstat_hits$gene_symbol[i] <- parts[1]
#   }
# }

f_test_data <- left_join(f_test_data, fstat_hits, by = "id")
f_test_data$highlight <- ifelse(f_test_data$id %in% genesymbol_hits, "Hit", "Not Hit")
f_test_data$label <- ""
#f_test_data[f_test_data$id == "BTK_1",]$label <- "BTK_1"
#f_test_data[f_test_data$id == "BTK_2",]$label <- "BTK_2"
f_test_data[f_test_data$id == "WASHC2C_1",]$label <- "WASHC2C_1"
f_test_data[f_test_data$id == "COG8_1",]$label <- "COG8_1"
#f_test_data[f_test_data$id == "BRAF_1",]$label <- "BRAF_1"
f_test_data[f_test_data$id == "MTOR_3",]$label <- "MTOR_3"
f_test_data[f_test_data$id == "FYN_1",]$label <- "FYN_1"
f_test_data[f_test_data$id == "ARID1A_1",]$label <- "ARID1A_1"
f_test_data[f_test_data$id == "VPS18_2",]$label <- "VPS18_2"
f_test_data[f_test_data$id == "SMARCD1_1",]$label <- "SMARCD1_1"
f_test_data[f_test_data$id == "SCYL1_1",]$label <- "SCYL1_1"
f_test_data[f_test_data$id == "CSK_1",]$label <- "CSK_1"


color_palette_unconventional <- c(
  "WASHC2C_1"= "#C79EFA",  
  "COG8_1"= "#FFB7B2", 
  "MTOR_3"= "#FFD6A5", 
  "FYN_1"=  "#B5EAD7", # Soft Mint
  "ARID1A_1"= "#CAFFBF", 
  "VPS18_2"= "#FAE1DD", 
  "SMARCD1_1"= "#B28DFF",
  "SCYL1_1"= "#A0C4FF",
  "CSK_1"= "#FDFFB6"
)



library(ggrepel)
pdf(file = "fstat_scatter_null_melt.pdf", width = 6.31, height = 5.27)
ggplot(f_test_data, aes(x = fStat, y = f_statistic, color = highlight, size = highlight, alpha = highlight)) +
  geom_point() +
  scale_color_manual(values = c("Hit" = "grey30", "Not Hit" = "grey80")) +
  scale_size_manual(values = c("Hit" = 3, "Not Hit" = 1)) +
  scale_alpha_manual(values = c("Hit" = 1, "Not Hit" = 0.5)) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "#1F78B4") + #padj of 0.05
  theme_minimal() +
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  guides(size = "none", alpha = "none", color = "none")+
  geom_label_repel(aes(label = label),max.overlaps = 660000, min.segment.length = 3) +
  geom_point(data = subset(f_test_data, label != ""), color = "#1F78B4", size = 5) +
  geom_point(data = subset(f_test_data, label != ""), color = "white", size = 4.9) +
  geom_point(data = subset(f_test_data, label != ""), color = color_palette_unconventional, size = 3) +
  labs(title = "F-Statistic Comparison",
       x = "F-Statistic, melt difference",
       y = "F-Statistic, abundance comparison",
       color = "Gene Symbol Status",
       size = "Gene Symbol Status",
       alpha = "Gene Symbol Status") +
  theme_pubclean()
dev.off()

f_test_data_significant <- f_test_data[#f_test_data$id %in% genesymbol_hits &
                                         f_test_data$adjusted_p_values < 0.05 ,] %>%
  select(., starts_with(c("gene_sym", "f_stat", "p_val", "adjusted_", "id")))

write.csv(f_test_data_significant, "f_test_data_significant_null.csv")

