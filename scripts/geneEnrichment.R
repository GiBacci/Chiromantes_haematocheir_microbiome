# ---- loadGeneData ---------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggalluvial)
  library(ggrepel)
  library(DBI)
  library(ggVennDiagram)
  library(GO.db)
})
source("scripts/utils.R")

# Font family
f.family <- "sans"

geneCount <- "data/gene_count.rds"
geneCount <- readRDS(geneCount)

# load EC2GO
ec2go <- read.table("data/ec2go.txt", sep = ";",
                    strip.white = T, skip = 2,
                    col.names = c("EC", "GO"),
                    quote = "") %>%
  separate(EC, into = c("EC", "description"), sep = " > ")


# ---- hTestGene ---------------------- 
# Hypergeometric test
go.enrich <- geneCount %>%
  mutate(
    p.value = case_when(
      p_obs > p_exp ~ phyper(q = x - 1, m = m, n = n, k = k, lower.tail = F),
      p_obs < p_exp ~ phyper(q = x, m = m, n = n, k = k, lower.tail = T),
      T ~ NA_real_
    ),
    adj.p.value = p.adjust(p.value, "BH")
  )

# Get significant GO terms
#  p.value < 0.05
#  |FC| > 1
sign.go <- go.enrich %>%
  filter(adj.p.value < 0.05,
         abs(log2.fold.enrichment) >= 1)


# Building table S9
table.s9 <- sign.go %>%
  mutate(Term = Term(GO),
         Reaction = Definition(GO)) %>%
  arrange(GO)


# ---- plotGeneEnrich --------------------
# Function for mixing colors 
# base on RGB value
mixCol <- function(col){
  if(any(is.na(col)))
    return(NA)
  
  if(length(col) == 1)
    return(col)
  
  x <- rowMeans(col2rgb(col))
  res <- rgb(x[1], x[2], x[3], 
             maxColorValue=255)
  paste0(res, "FF")
}

# Defining base colors (the same used 
# throughout the paper)
col <- ggsci::pal_npg()(10)[1:3]
names(col) <- LETTERS[1:3]

# Splitting dataset based on groups
sign.go.split <- sign.go %>%
  left_join(ec2go, by = "GO") %>%
  mutate(descrfc = paste0(description, sign(log2.fold.enrichment))) %>%
  split(.$Cluster)
# Setting names for plotting
names(sign.go.split) <- paste("Cluster", names(sign.go.split))

# Building Venn diagram
venn <- sign.go.split %>%
  map("descrfc") %>%
ggVennDiagram(color = "white", label = "count")

# Detecting geoms
venn.geoms <- map(venn$layers, "geom") %>%
  map_chr(function(x) class(x)[1])

poly.layer <- which(venn.geoms == "GeomPolygon")
text.layer <- which(venn.geoms == "GeomText")
label.layer <- which(venn.geoms == "GeomLabel")

# Assign mixed colors
col <- venn$layers[[poly.layer]]$data$group %>%
  unique() %>%
  set_names() %>%
  strsplit("") %>%
  map(function(x) {
    res <- mixCol(col[x])
    setNames(res, NULL)
    }) %>%
  unlist()

# Venn diagram customizations
venn$layers[[text.layer]]$aes_params$size <- 2.8
venn$layers[[text.layer]]$aes_params$fontface <- "plain"

venn$layers[[label.layer]]$aes_params$alpha <- 0
venn$layers[[label.layer]]$aes_params$size <- 2.8

venn$layers[[poly.layer]]$mapping <- aes(group=group, fill=group)
# venn$layers[[poly.layer]]$data <- venn$layers[[2]]$data %>%
#   mutate(group = fct_relevel(group, "A", "B", "C"))

# Final plot
venn <- venn + 
  scale_fill_manual(values = col) +
  myTheme +
  theme_void(base_size = 8, base_family = f.family) +
  theme(legend.position = "none")

# Size of clusters (number of unique GO terms)
bars <- sign.go.split %>%
  map("descrfc") %>%
  map(unique) %>%
  map_int(length) %>%
  enframe("group", "size") %>%
ggplot(aes(x = size, y = group)) +
  geom_col(aes(fill = group)) +
  theme_classic(base_size = 8, base_family = f.family) +
  myTheme +
  scale_x_continuous(expand = expansion(mult = c(0,0.1))) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank()) +
  ggsci::scale_fill_npg() +
  xlab("N° of GO terms")


# Get venn groups
venn.grp <- venn$layers[[poly.layer]]$data %>%
  dplyr::select(group, count) %>%
  unique()

# Top GO (based on enrichment fold change)
top_go <- bind_rows(sign.go.split) %>%
  # group by term
  group_by(descrfc) %>%
  # get set/s in which a term is present
  mutate(set = paste(unique(Cluster), collapse = ",")) %>%
  group_by(set, .add = T) %>%
  # compute mean FC for each set 
  do(mean_cl_boot(.$log2.fold.enrichment)) %>%
  group_by(set) %>%
  # count the number of terms in each set
  # and match it with the groups in the 
  # venn diagram
  mutate(count = n()) %>%
  left_join(venn.grp, by = "count") %>%
  # Fine tuning text
  mutate_at("descrfc", str_remove, "(activity)?[0-9]$") %>%
  mutate_at("descrfc", str_remove, "^GO:") %>%
  # select top N terms
  top_n(10, abs(y)) %>%
  ungroup() %>%
  arrange(abs(y)) %>%
  # Fine tuning factor order
  mutate_at("descrfc", fct_inorder) %>%
  mutate_at("group", fct_relevel, "ABC", after = Inf) %>%
ggplot(aes(x = descrfc, y = y)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), size = .3,
                width = .5) +
  geom_col(fill = "steelblue") +
  facet_grid(group ~ ., scales = "free") +
  theme_bw(base_size = 8, base_family = f.family,
                base_line_size = .25) +
  coord_flip() +
  scale_y_continuous(expand = expansion(add = c(0, .3))) +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank()) +
  ylab("Fold-change enrichment (log2)")

# Use venn diagram groups as facet label
venns <- unserialize(serialize(venn, NULL))
venns$layers <- venns$layers[poly.layer]
venns$layers[[1]]$data <- top_go$data$group %>% 
  map_df(function(l){
    venns$layers[[1]]$data %>%
      mutate(group.col = ifelse(group == l, as.character(group), NA),
             panel = l)
  })
venns$layers[[1]]$mapping <- aes(group = group, fill = group.col)

# Multipanel venns coloring only the selected
# group (we could re-use the plot for the heatmap)
venns <- venns + 
  scale_fill_manual(values = col, na.value="gray90") +
  facet_wrap(~ panel, ncol = 1) +
  theme(strip.background = element_blank(),
        strip.text = element_blank())

# Building final plot
top_venn <- top_go + 
  (venns + theme(panel.spacing.y = unit(2, "lines"))) +
  plot_layout(widths = c(2,1))


##### Plot enriched genes #####
genes <- readRDS("data/gene_asv_matrix.rds")
group <- readRDS("output/groups.rds") %>%
  filter(!is.na(group))

# Joining GO terms significantly enriched:
#  to avoid overlaps, the direction of the fold changes
#  (-1 or +1) is added to the go Terms and the log fold change
#  are averaged.
#  For each ASV the median number of copies is computed in
#  each GO term.
collapsed <- genes %>%
  filter(GO %in% sign.go$GO) %>%
  # Joining with significant table
  left_join(dplyr::select(sign.go, GO, group.go = Cluster, lfc = log2.fold.enrichment), 
            by = "GO") %>%
  # Adding fold change direction
  mutate(descr = paste(description, sign(lfc), sep = "_")) %>%
  # Group by ASV and GO term (with FC)
  group_by(taxon, descr) %>%
  # Averaging fold change and median of copy number
  summarise(group.go = paste(sort(unique(group.go)), collapse = ","),
            mfc = mean(lfc),
            median_taxon_ab = median(taxon_function_abun) * sign(mfc))

# Adding ASV clusters and removing those with
# no cluster attribution (cluster 4 has no 16S)
collapsed <- collapsed %>%
  left_join(dplyr::select(group, taxon = asv, 
                          group.taxon = group), by = "taxon") %>%
  filter(!is.na(group.taxon))

# Relabel function for matching GO terms
# groups with the venn groups computed by 
# the 'ggVennDiagram' function.
relabel <- function(x){
  i <- strsplit(x, ",")
  sapply(i, function(y){
    paste(LETTERS[as.numeric(y)], collapse = "")
  })
}

# Re-lable groups and ordering levels to match the
# venn plot
collapsed <- collapsed %>%
  ungroup() %>%
  mutate(descr = fct_reorder(descr, mfc, .fun = mean, .desc = T),
         group.go = fct_relabel(group.go, relabel),
         group.go = fct_relevel(group.go, levels(venns$layers[[1]]$data$panel)))

# HeatMap of copy_numbers after 
# hierarchical clustering based on kendall
# correlation
mat <- collapsed %>%
  mutate(median_taxon_ab = abs(median_taxon_ab)) %>%
  dplyr::select(taxon, descr, median_taxon_ab) %>%
  spread(descr, median_taxon_ab, fill = 0) %>%
  data.frame(row.names = 1) %>%
  t() %>% scale() %>%
  cor(method = "kendall")
mat <- as.dist((1 - mat)^2)
ord.asv <- hclust(mat, "average") %>%
  as.dendrogram() %>%
  order.dendrogram()

# Color limits (simmetric around zero)
fill.lim <- log10(max(abs(collapsed$median_taxon_ab))) * c(-1,1)
# heatmap
heatmap_copies <- collapsed %>%
  mutate(median_taxon_ab = log10(abs(median_taxon_ab)) * sign(median_taxon_ab),
         taxon = factor(taxon, levels = labels(mat)[ord.asv])) %>%
  ggplot(aes(x = taxon, y = descr, fill = median_taxon_ab)) +
  geom_tile(color = NA) +
  facet_grid(group.go ~ as.factor(group.taxon),
             scales = "free", space = "free") +
  theme_minimal(base_size = 8, base_family = f.family,
                base_line_size = .25) +
  scale_fill_distiller(palette = "RdBu", limit = fill.lim,
                       labels = function(x) abs(as.numeric(x))) +
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = "transparent"),
        legend.position = "left",
        # legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(fill = guide_colorbar(label.position = "left",
                               barwidth = unit(.5, "lines"),
                               barheight = unit(15, "lines"),
                               frame.colour = "black",
                               frame.linewidth = .25,
                               ticks.colour = "black",
                               title.position = "left",
                               title.theme = element_text(angle = 90, size = 8),
                               title.hjust = .5,
                               title = bquote(median~taxon~abundance~(log[10])))) +
  xlab("Enriched ASVs") +
  ylab("Depleted | Enriched GO terms")

# Table S8
count.terms <- collapsed %>%
  group_by(group.go, descr) %>%
  summarise(mfc = sign(mean(mfc))) %>%
  group_by(mfc, .add = T) %>%
  dplyr::count()

mat <- matrix(nrow = nrow(count.terms), ncol = 3)
for(i in 1:3){
  l <- LETTERS[i]
  rows <- grep(l, as.character(count.terms$group.go))
  mat[rows, i] <- "X"
}
mat[is.na(mat)] <- ""
mat <- data.frame(mat, count.terms[,-1])
colnames(mat) <- c(paste0("C", 1:3), "Effect", "N")

mat <- mat %>%
  mutate(Effect = case_when(
    Effect == -1 ~ "Depleted",
    Effect == 1 ~ "Enriched"
  )) %>%
  pivot_wider(names_from = Effect, 
              values_from = N, 
              values_fill = 0)

totals <- t(sapply(1:3, function(i){
  tot <- colSums(mat[mat[[i]] == "X", 4:5])
  labs <- rep("", 3)
  labs[i] <- "Tot"
  c(labs, tot)
}))
totals <- data.frame(totals)
names(totals) <- names(mat)
for(i in 4:5){
  totals[[i]] <- as.numeric(totals[[i]])
}
table.s8 <- rbind(mat, totals)

# Fold change bars
fc_bars <- collapsed %>%
  group_by(descr, group.go) %>%
  summarise(mfc = mean(mfc)) %>%
  ggplot(aes(y = descr, yend = descr, x = 0, xend = mfc)) +
  geom_segment(color = "steelblue") +
  geom_vline(xintercept = 0, size = .25, color = "black") +
  facet_grid(group.go ~ .,
             scales = "free", space = "free") +
  theme_minimal(base_size = 8, base_family = f.family,
                base_line_size = .25) +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.x = element_line(),
        axis.title.y = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = "transparent"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  xlab("Fold-change enrichment (log2)")

# Copy venn diagrams to match only cluster 1, 2, and 3
venns.main <- unserialize(serialize(venns, NULL))
venns.main$layers[[1]]$data <- venns.main$layers[[1]]$data %>%
  filter(panel %in% c("A", "B", "C"))
venns.main <- venns.main + facet_wrap(~ panel, nrow = 1) 


# Final format (it needs a manual fine tune)
layout <- "
EEEEEEEEFFF
EEEEEEEEFFF
AAAAAAAAFFF
BBBBBBBBCCD
BBBBBBBBCCD
BBBBBBBBCCD
BBBBBBBBCCD
BBBBBBBBCCD
BBBBBBBBCCD
BBBBBBBBCCD
BBBBBBBBCCD
BBBBBBBBCCD
BBBBBBBBCCD
BBBBBBBBCCD
BBBBBBBBCCD
BBBBBBBBCCD
"
figure.s8 <- (venns.main + theme(panel.spacing.x = unit(4, "lines"))) + 
  (heatmap_copies + labs(tag = "c")) + 
  (fc_bars + labs(tag = "d")) + 
  (venns + theme(panel.spacing.y = unit(1, "lines"))) +
  (bars + labs(tag = "a")) + 
  (venn + coord_equal(clip = "off") + labs(tag = "b")) +
  plot_layout(design = layout) &
  theme(plot.tag = element_text(vjust = -3, face = "bold", 
                                size = myTheme$text$size + 1))
