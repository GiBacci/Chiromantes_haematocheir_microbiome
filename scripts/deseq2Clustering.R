# ---- loadDataClust ---------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(patchwork)
  library(ggbeeswarm)
  library(phyloseq)
  library(cluster)
  library(ggtree)
  library(ggsci)
  library(data.tree)
  library(DBI)
})
source("scripts/utils.R")
source("scripts/hypergeometricTest.R")
source("scripts/sunburstPlot.R")

# Reading experiment
exp <- readRDS("output/exp.rds")
# Removing singletons
exp <- filter_taxa(exp, function(x) sum(x) > 1, T)

# get TAX table
tax <- tax_table(exp) %>%
  as("matrix") %>%
  data.frame() %>%
  rownames_to_column("asv")

# ---- modelDefinition -------------------------

# Fitting models and save
dds.out <- "output/deseq2_fits.rds"
# if output exists load it otherwise
# fit models from scratch
if(!file.exists(dds.out)){
  dds <- phyloseq_to_deseq2(exp, ~ env + site)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  
  fits <- c("parametric", "local", "mean") %>%
    set_names() %>%
    map(estimateDispersions, object = dds)
  
  saveRDS(fits, dds.out)
}else{
  fits <- readRDS(dds.out)
}

# Plot fit estimates
disp.est <- fits %>% 
  map(mcols) %>%
  map_df(as_tibble, rownames = "asv", .id = "fitType")

# Compute median absolute log residual
# between estimated dispersion and
# fitted dispersion 
labs <- disp.est %>%
  filter(!is.na(dispGeneEst)) %>%
  group_by(fitType) %>%
  summarise(mar = median(abs(log(dispGeneEst) - log(dispFit))),
            x = max(baseMean),
            y = min(dispGeneEst)) %>%
  ungroup()

# Plot dispersions
disp.all <- disp.est %>%
  dplyr::select(baseMean, dispGeneEst, dispersion, dispFit, dispOutlier, fitType) %>%
  pivot_longer(c(-baseMean, -fitType, -dispOutlier),  names_to = "type") %>%
  mutate(type = fct_recode(type, 
                           final = "dispersion",
                           estimated = "dispGeneEst",
                           fitted = "dispFit"),
         outlier = dispOutlier & type == "final")
# Get outliers
disp.out <- filter(disp.all, outlier)
# Get fitting
disp.fitted <- filter(disp.all, type == "fitted")

# Plot dispersions
figure.s6 <- disp.all %>%
  filter(!outlier, type != "fitted") %>%
  ggplot(aes(x = baseMean, y = value, 
             color = type, size = type)) +
  geom_point(alpha = .5) +
  # Add outliers
  geom_point(data = disp.out, shape = 1, size = .8, alpha = .5) +
  # Add fitted model
  geom_line(data = disp.fitted, show.legend = F) +
  # Adding labels (meadin absolute residual)
  geom_text(data = labs, aes(x = x, y = y + 1,
                             label = round(mar, digits = 3)),
            inherit.aes = F, hjust = 1, vjust = 0, size = 3) +
  scale_x_log10() +
  scale_y_log10(expand = c(0.1, 0.1)) +
  # Plot grid
  facet_wrap(~ fitType, ncol = 1, strip.position = "right",
             labeller = labeller(fitType = function(x) paste("Fit type:", x))) +
  # Adjusting scales manually
  scale_color_manual(values = c("gray10", "dodgerblue3", "firebrick3")) +
  scale_size_manual(values = c(.2, .4, .5)) +
  # Suppress legend for size and shape
  guides(color = guide_legend(title = NULL,
                              override.aes = list(size = 2, shape = 19)),
         size = FALSE, shape = FALSE) +
  theme_bw(base_size = 8, base_family = "sans") +
  theme_publication() +
  xlab("Mean of normalized counts") +
  ylab("Dispersion") +
  coord_cartesian(ylim = c(1, 200))

# Local fit is the best fit (the one with the
# lower meadin absolute residual)
dds <- fits$local

# ---- LRTtest -------------------------
# LRT test. If output file exists do
# not test again
out <- "output/deseq2_lrt.rds"
if(!file.exists(out)){
  # LRT test
  envs <- nbinomLRT(dds, reduced = ~ site)
  saveRDS(envs, out)
}else{
  envs <- readRDS(out)
}
table.s7 <- results(envs, alpha = 0.05) %>%
  as("data.frame") %>%
  rownames_to_column("asv")

# Significant table with taxonomic annotations
sign.asv <- results(envs, alpha = 0.05) %>%
  as("data.frame") %>%
  rownames_to_column("asv") %>%
  filter(!is.na(padj), padj <= 0.05) %>%
  left_join(tax, by = "asv")

# Pulling significant ASVs
asvs <- sign.asv %>%
  pull(asv)

asvs.count <- prune_taxa(asvs, exp) %>%
  tax_table() %>%
  as("matrix") %>%
  as.data.frame() %>%
  pull("domain") %>%
  table()

# Summary statistics
sign.ab <- prune_taxa(asvs, exp) %>%
  sample_sums()
all.ab <- sample_sums(exp)
sign.mean <- (sign.ab / all.ab) * 100
sign.se <- sd(sign.mean) / sqrt(length(sign.mean))
sign.mean <- mean(sign.mean)
sign.ab <- (sum(sign.ab) / sum(all.ab)) * 100

taxa.percent <- (length(asvs) / ntaxa(exp)) * 100

# ---- deseqclustering --------------------

# Getting matrix for clustering
mat <- vst(envs, fitType='local') %>%
  assay()
meta <- colData(envs) %>%
  as("data.frame")


# Getting mean by groups
grp <- meta$env
means <- apply(mat[asvs,], 1, function(x){
  by(x, grp, mean)
})
grp <- rownames(means)
means <- apply(means, 2, scale)
rownames(means) <- grp

# Computing distances based on kendall
# correlation
d <- as.dist((1 - cor(means, method = "kendall"))^2)

### Diana clustering ###########################
cl <- diana(d, diss = T, stand = F)

# Grouping based on divisive coefficient
hc <- as.hclust(cl)
group <- cutree(hc, h = cl$dc)
# Removing clusters with less than 15 ASVs
group[group %in% names(table(group))[table(group) <= 15]] <- NA
# Building data frame
group <- tibble(asv = names(group),
                group = dense_rank(group)) %>%
  left_join(tax, by = "asv")

saveRDS(group, "output/groups.rds")
# Plot clustering ---------------------------------------------------------

# Getting clades based on most recent 
# common ancestor algorithm (MRCA)
g <- group %>%
  split(.$group) %>%
  map(pull, asv)

# Building tree and grouping clades
# based on clustering
p <- ggtree(hc, size = .6, lineend = "square")
clades <- sapply(g, function(n) MRCA(p, n))


p <- groupClade(p, clades, group_name='subtree')
# MOdifying labels
p$data <- p$data %>%
  mutate(subtree = fct_relabel(subtree, ~paste("Cluster", .x))) %>%
  mutate(subtree = fct_recode(subtree, Removed = "Cluster 0"))

# Definig color for clustering
groups <- p$data %>%
  dplyr::count(subtree)
col <- c("black", pal_npg()(10)[1:(nrow(groups)) - 1])
names(col) <- pull(groups, subtree)

# Plotting with heatmap
p$data$x <- ifelse(p$data$x > 0, 0, p$data$x)
p <- p + aes(color=subtree) +
  geom_treescale(x=min(p$data$x) + .1, 
                 y=max(p$data$y) - 10, 
                 offset = 2,
                 family = "sans", fontsize = 3,
                 width = .2)
figure.5 <- gheatmap(p %<+% group, t(means), color = NA,
                       colnames_offset_y = -4,
                       offset = .04,
                       width = .5,
                       font.size = 2,
                       family = "sans",
                       colnames_angle = 90,
                       hjust = 1) +
  geom_tippoint(aes(shape = domain, x = x+0.05), size = 1, 
                stroke = 0, show.legend = F) +
  scale_color_manual(values = col) +
  scale_shape_manual(values = c(NA, 19)) +
  scale_fill_distiller(palette = "RdBu", limits = c(-2.5, 2.5)) +
  scale_y_continuous(expand = expansion(add=c(10,0))) +
  scale_x_continuous(expand = c(0,0)) +
  theme_minimal(base_size = 8, base_family = "sans") +
  theme_publication(grid = F, axis.line = F,
                    legend.position = "right",
                    legend.direction = "vertical",
                    ticks.len = 0) +
  theme(plot.margin = unit(c(.2,.2,.2,.2), "lines"),
        axis.text = element_blank()) +
  guides(color = guide_legend(title.position = "top",
                              title = "Cluster",
                              override.aes = list(size = 2)),
         fill = guide_colorbar(title = "Scaled mean\nabundance (VST)",
                               barheight = unit(10, "lines"),
                               barwidth = unit(.5, "lines")))


# ---- pairwiseWilc ----------------------------------------------

# Function for computing pairwise wilcox test
# and returning letters for significant contrasts
#    data : data frame
#    value : name of the column with values to be tested
#    grp : name of the column with grouping factor
compareGroups <- function(data, value, grp){
  library(multcompView)
  
  # Getting formula and data
  frm <- as.formula(paste(value, grp, sep = "~"))
  
  x <- get(value, data)
  g <- get(grp, data)
  
  # testing
  test <- pairwise.wilcox.test(x, g, p.adjust.method = "BH")
  
  # Getting p-value matrix
  pvals <- data.frame(test$p.value) %>%
    rownames_to_column(var = "A")
  
  # Formatting data frame with contrasts
  # to be used in "multcompLetters2" function
  test <- pvals %>%
    gather("B", "p", -A) %>%
    filter(!is.na(p)) %>%
    mutate(name = paste(A,B, sep = "-")) %>%
    dplyr::select(name, p) %>%
    deframe()
  
  # Setting contrasts evaluation
  by.col <- setNames(c("A"), grp)
  test <- multcompLetters2(frm, test, data)
  
  # Joining letters and p-values
  test <- enframe(test$Letters, 
                  name = "env", value = "label") %>%
    left_join(pvals, by = by.col)
  
  # Computinf y-values for plotting letters
  res <- by(x, g, max, simplify = T)
  
  # Retunring a data frame
  enframe(res, name = grp, value = value) %>%
    left_join(test, by = grp)
}

# Getting all significant contrasts
# (supplementary materials)
samp <- meta %>%
  group_by(env, type) %>%
  summarise()

# Getting data divided by clusters
data <- g %>%
  map_df(function(i){
    means[,i] %>%
      as_tibble(rownames = "env") %>%
      gather("asv", "value", -env)
  }, .id = "cluster") %>%
  group_by(cluster) %>%
  mutate(n = length(unique(asv)),
         cluster.lab = paste0("Cluster ", cluster, " (", n, " ASVs)")) %>%
  ungroup() %>%
  left_join(samp, by = "env") %>%
  mutate_at("env", fct_inorder) %>%
  mutate_at("type", fct_inorder)

# Testing clusters with "compareGroups"
# funciton
test <- data %>%
  data.frame() %>%
  split(.$cluster.lab) %>%
  map_df(compareGroups, value = "value", grp = "env",
         .id = "cluster.lab") %>%
  left_join(samp, by = "env") %>%
  mutate_at("type", fct_inorder)


# Plotting p-values and boxplots
# Pretty complex
a <- ggplot(data, aes(x = env, y = value)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1, 
               aes(fill = cluster), size = .23) +
  geom_text(data = test, aes(label = label),
            vjust = -1, size = 3) +
  facet_grid(cluster.lab ~ type, space = "free_x", scales = "free_x",
             labeller = labeller(type = function(x){
               case_when(
                 x == "biotic" ~ "Crab's organs",
                 T ~ "Environmental samples"
               )
             })) +
  scale_y_continuous(expand = expansion(add = c(0, .8))) +
  scale_fill_npg() +
  myTheme +
  ylab("Mean scaled abundance (VST)") +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        panel.grid.major.x = element_blank())


ord <- data %>% pull(env) %>% levels()
b <- test %>%
  dplyr::select(-value, -label, -type) %>%
  gather("env.2", "p", -cluster.lab, -env) %>%
  filter(!is.na(p)) %>%
  mutate(sign = case_when(
    p < 0.01 ~ "p < 0.01",
    p < 0.05 ~ "p < 0.05",
    TRUE ~ "p >= 0.05"
  )) %>%
  mutate_at(c("env", "env.2"), fct_relevel, ord) %>%
  ggplot(aes(x = env, y = env.2, fill = sign)) +
  geom_tile(color = "white") +
  facet_wrap(~ cluster.lab, ncol = 1) +
  coord_cartesian(expand = F) +
  myTheme +
  scale_y_discrete(position = "right", expand = expansion(mult = .1)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.direction = "vertical",
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.key.height = unit(.4, "lines"),
        legend.key.width = unit(.7, "lines")) +
  scale_fill_manual(values = c("firebrick", "royalblue", "gray"))

# Final Plot
figure.s7 <- cowplot::plot_grid(a, b + theme(plot.margin = unit(c(1.2, .2,.4,.2), "lines")), 
                                labels = c("a", "b"), rel_widths = c(3.8, 2),
                                label_size = myTheme$text$size + 1, axis = "bt")


# ---- taxaEnrich -----------------------------------------

# A plot similar to panle a of the previous 
# figure but suitable for the main text
# togheter with enrichment analysis

# Cluster 1 : most abundant in MG and HG.
#             A double asterisk is reported
#             on top of the two groups.
#
# Cluster 2 : most abundant in GI. A double
#             asterisk will be reported on
#             top of it
#
# Cluster 3 and 4 : environmental clusters.
#                   The first in most abundant in W
#                   and WD (water cluster), whether 
#                   the second in L and S (earth cluster).
lvl <- data %>% pull(cluster.lab) %>% unique()

# Adding asterisks only for contrasts described above
# (other contrasts are reported in the previous figure)
asterisks <- data.frame(
  cluster.lab = c(lvl[1], lvl[1], lvl[2], lvl[3], lvl[3], lvl[4], lvl[4]),
  env = c("MG", "HG", "GI", "W", "WD", "L", "S"),
  x = c(1, 2, 4, 3, 4, 1, 2),
  label = c("**", "**", "**", "**", "**", "**", "**"),
  type = c("biotic", "biotic", "biotic", "abiotic", "abiotic", "abiotic", "abiotic")
)

# Joining with data for y evaluation
lab <- data %>%
  group_by(env, cluster.lab) %>%
  summarise(value = mean(value)) %>%
  right_join(asterisks, by = c("env", "cluster.lab"))

# Adding contrasts between most represented
# groups
asterisks <- data.frame(
  env = c("MG", "W"),
  cluster.lab = c(lvl[1], lvl[3]),
  value = 0,
  x = c(1.5, 3.5),
  label = c("**", "*"),
  type = c("biotic", "abiotic")
)


# Plotting
a <- ggplot(data, aes(x = env, y = value, fill = cluster.lab)) +
  stat_summary(geom = "col", fun = mean, show.legend = F) +
  geom_text(data = lab, aes(x = x, label = label), vjust = 0) +
  geom_text(data = asterisks, aes(x = x, label = "___"), vjust = 0.3, size = 3) +
  geom_text(data = asterisks, aes(x = x, label = label), vjust = 1.5) +
  facet_grid(cluster.lab ~ type, space = "free_x", scales = "free_x",
             labeller = labeller(type = function(x){
               case_when(
                 x == "abiotic" ~ "Environmental\nsamples",
                 x == "biotic" ~ "Crab's\norgans"
               )
             }, 
             cluster.lab = function(x){
               str_replace(x, " \\(", "\n(")
             })) +
  geom_hline(yintercept = 0, size = .25) +
  scale_y_continuous(expand = expansion(add = c(.1, .5))) +
  myTheme +
  scale_fill_npg() +
  theme(panel.spacing.x = unit(.1, "lines"),
        # legend.position = "none",
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  ylab("Mean scaled abundace (VST)")


# Enrichment
clusters <- data.frame(asv = taxa_names(exp)) %>%
  left_join(group, by = "asv") %>%
  pull(group)

cat <- buildDataSet(mat, fct = meta$env, grp = clusters)
t.table <- tax_table(exp) %>%
  as("matrix") %>%
  as_tibble()

# Hierarchical enrichment
enrich <- as.list(t.table)[-1] %>%
  imap_dfr(function(ctg, name){
    res <- hpostSingle(cat, as.factor(ctg)) %>%
      filter(adj.p.value < 0.05,
             log2.fold.enrichment > 1)
    
    tibble(
      res, 
      t.table[match(res$categories, ctg),] %>%
        mutate_at(vars(name:species, -name), 
                  function(x) rep(NA, length(x)))
    )
    
  }, .id = "tax.lvl")

table.s8 <- enrich %>% 
  dplyr::rename(
    "Level" = "tax.lvl",
    "Cluster" = "grp",
    "Taxon" = "categories",
    "x" = "q",
    "p_obs" = "fraction.obs",
    "p_exp" = "fraction.pop",
    "Amplicon" = "superdom"
  ) %>%
  dplyr::select(-exp.q) %>%
  dplyr::select(Amplicon, Level, Taxon, Cluster, x, k, n, m, p_obs, 
         p_exp, log2.fold.enrichment, p.value, adj.p.value)

# ---- enrichmentPlots --------------------------
data <- tax %>%
  filter(asv %in% asvs) %>%
  dplyr::select(-asv)

# Removing taxa without significant enrichment
data <- enrich %>%
  dplyr::select(tax.lvl, categories) %>%
  unique() %>%
  apply(1, function(d){
    col.i <- which(names(data) == d[1])
    row.i <- data[,col.i] == d[2]
    res <- data[row.i, 1:col.i]
    unique(res)
  }) %>%
  bind_rows() %>%
  unique()
rownames(data) <- NULL
data <- data[,-1]
extr <- apply(data, 1, function(x) !all(is.na(x)))
data <- data[extr,]

# ensuring all taxa were correctly found
diffs <- setdiff(enrich$categories, na.omit(as.vector(as.matrix(data))))
if(length(diffs) != 0){
  print("[WARNING] Enriched taxa were not correctly detected")
}

# Building sunburst plot
res.raw <- getSunburstData(data, collapse = T)

# Defining cluster colors
col <- pal_npg()(10)[1:4]

# Function for mixing colors
mixCol <- function(col){
  if(any(is.na(col)))
    return(NA)
  
  if(length(col) == 1)
    return(col)
  
  x <- rowMeans(col2rgb(col))
  res <- rgb(x[1], x[2], x[3], maxColorValue=255)
  paste0(res, "FF")
}

correctLabels <- function(data, label){
  data %>%
    mutate(across({{ label }}, str_remove, "Allorhizobium-Neorhizobium-Pararhizobium-")) %>%
    mutate(across({{ label }}, str_remove, "Cryptococcus sp XVII ")) %>%
    mutate(across({{ label }}, str_replace, "Incertae Sedis", "i.s.")) %>%
    mutate(across({{ label }}, str_replace, "_Incertae sedis", " i.s."))
}

# Adding clusters
res <- res.raw %>%
  left_join(enrich %>% dplyr::select(grp, categories), 
            by = c("label" = "categories")) 

# Interpolate colors if a taxa was
# enriched in two or more clusters
res <- res %>%
  mutate(col = col[grp]) %>%
  group_by(label) %>%
  mutate(col = mixCol(unique(col)),
         lab = paste(unique(grp), collapse = "&")) %>%
  ungroup() %>%
  unique()

# Get color vector
col <- na.omit(deframe(res %>% dplyr::select(lab, col)))
col.names <- unique(names(col))
col <- setNames(unique(col), col.names)

# Get labels with significant enrichment
lab <- res %>%
  filter(!is.na(grp)) %>%
  pull(label) %>%
  unique()

# Get external labels
lab1 <- res.raw %>% 
  filter(last == T,
         label %in% lab) %>%
  correctLabels(label)

# Get row labels
lab2 <- res.raw %>%
  filter(xmin == 1,
         label %in% lab) %>%
  correctLabels(label)

# Get labels big enough to be displayed
bigLabs <- res.raw %>%
  mutate(char.len = nchar(label)) %>%
  filter(ymax - ymin >= char.len - 1,
         xmin > 1,
         last == F,
         label %in% lab) %>%
  correctLabels(label)

lab3 <- res.raw %>%
  filter(xmin > 1,
         last == F,
         label %in% lab,
         !label %in% pull(bigLabs, label)) %>%
  correctLabels(label)

# Final version of plot
text.size <- 1.8
ext.col <- "gray50"
key_glyph <- get_draw_labels(8)
b <- res %>% mutate(lab = ifelse(lab == "NA", "n.s.", lab),
                    lab = fct_relevel(lab, "1", "2", "3", "4", "1&2", "2&4", "3&4", "n.s.")) %>%
ggplot(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_rect(aes(fill = lab), color = "white", show.legend = T) +
  polar_labels(lab1, alignment = "perpendicular", size = text.size, show.legend = F) +
  polar_labels(lab2, alignment = "tangent", position = "inside", 
               size = text.size, show.legend = F, use.shadowtext = T,
               color = "black", bg.colour = "white",
               fontface = "bold") +
  polar_labels(bigLabs, alignment = "tangent", position = "inside", 
               size = text.size, show.legend = F, use.shadowtext = T,
               color = "black", bg.colour = "white",
               fontface = "bold") +
  polar_labels(lab3, alignment = "normal", position = "inside", 
               condensed = T, condensed.legend = T, size = text.size,
               key_glyph = key_glyph, use.shadowtext = T,
               color = "black", bg.colour = "white",
               fontface = "bold") +
  coord_polar("y", clip = "off") +
  xlim(-1, NA) +
  scale_fill_manual(values = c(col, n.s. = "gray90")) +
  theme_void(base_size = 8, base_family = "sans") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        # legend.title = element_blank(),
        legend.key.width = unit(.5, "lines"),
        legend.key.height = unit(.5, "lines"),
        plot.margin = unit(c(0, 0, 0, 0), "lines")) +
  guides(fill = guide_legend(title = "Clusters", ncol = 1, override.aes = list(color = NA),
                             title.position = "top"),
         label = guide_legend(title = NULL, nrow = 15, override.aes = list(fill = "transparent")))

# Assembling figure 6
figure.6 <- (a + theme(strip.text.x = element_text(angle = 90, hjust = 0))) + b +
  plot_layout(ncol = 2, widths = c(1, 4),
              guides = "collect") +
  plot_annotation(tag_levels = "a") & 
  theme(legend.position = 'bottom',
        plot.tag = element_text(vjust = -3, face = "bold", 
                                size = myTheme$text$size + 1))
