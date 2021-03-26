# ---- loadData ---------------------------------
suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(tidyverse)
  library(ggsci)
  library(patchwork)
  library(ggpubr)
  library(DESeq2)
  library(cowplot)
  library(iNEXT)
  library(UpSetR)
})
source("scripts/utils.R")

# Set parameters
#  binary : should alpha diversity be calculated using a binary index?
#  subsample : should sample be subsampled?
binary <- F
subsample <- NULL

# Function that separates the "data" column
# into amplicon and type
separateVar <- function(x){
  x %>% separate(data, into = c("amplicon", "type"), sep = "__")
}

# Normalising using DESeq2
exp <- readRDS("output/exp.rds")

# Percentage of taxa detected (presence/absence)
# Used to manually add results in the main text
res <- apply(tax_table(exp), 2, table)
res <- sapply(res, function(x) {
  perc <- round((x / sum(x)) * 100, 2)
  res <- paste0(x, " (", perc, "%)")
  setNames(res, names(x))
})

# Building DESeq2 object
dds <- phyloseq_to_deseq2(exp, ~ env + site)
dds <- estimateSizeFactors(dds, type = "poscounts")
# Normalise counts for sequencing depth:
#  Counts are divided by size factors obtained with
#  the "median ratio method".
otu_table(exp) <- counts(dds, norm = T) %>%
  t() %>% otu_table(taxa_are_rows = F)

# Subsampling if necessary (default to FALSE)
if(!is.null(subsample)){
  crabs <- sample(1:11, subsample)
  envs <- sample(1:5, subsample)
  exp <- subset_samples(exp, crab %in% c(crabs, 0))
  exp <- subset_samples(exp, sample %in% envs | type == "biotic")
}

# Format labels
labelFrmt <- function(x){
  case_when(
    x == "biotic" ~ "Crab's\norgans",
    x == "abiotic" ~ "Environmental\nsamples",
    grepl("\\|", x) ~ "Both",
    T ~ x
  )
}


# ---- alphaDiversity -----------------------------

# Launching interpolation/extrapolation only
# if output file dos not exist. To repeat
# this step (time consuming) simply remove the
# "iNEXT.out" file.
iNEXT.out <- "output/alpha.iNEXT.rds"
if(!file.exists(iNEXT.out)){
  # iNEXT package (interpolation / extrapolation)
  # Splitting experiment
  exp.rare <- splitBy(exp, 
                      sample = "env", 
                      taxa = "superdom", 
                      both = F, 
                      keep.zeroes = F)
  
  # Launching iNEXT with multiple threads:
  # 1. Setting up multithreading system based
  #    on operating system (windows: doSNOW, linux: doMC)
  threads <- 4
  if(Sys.info()['sysname'] == "Windows"){
    library(doSNOW)
    cl <- makeCluster(threads)
    registerDoSNOW(cl)
  }else{
    library(doMC)
    old <- getDoParWorkers()
    registerDoMC(threads)
  }
  
  # 2. Launching iNEXT
  rare <- foreach(e = exp.rare, .final = function(x) 
    setNames(x, names(exp.rare)),
    .packages = c("phyloseq", "iNEXT", "tidyverse")) %dopar% {
      # Get count matrix
      x <- otu_table(e) %>% as("matrix") %>% t() %>% round()
      # Get amplicon name
      amplicon <- tax_table(e)[,1] %>% unique() %>% as.vector()
      
      # Setting breaks according to amplicon counts:
      #   16S = 470000 max counts (300 breaks)
      #   ITS =  60000 max counts (100 breaks)
      breaks <- if(amplicon == "ITS"){
        seq(1, 60000, length.out = 100)
      }else{
        seq(1, 470000, length.out = 300) 
      }
      
      # Estimating rarefaction
      # A value of q equal to 2 corresponds to the inverse simpson
      # index.
      iNEXT(x, q = 2, datatype = "abundance", se = TRUE, 
            conf = 0.95, nboot = 30,
            size = breaks)
    }
  # 3. Closing mutithreading
  if(Sys.info()['sysname'] == "Windows"){
    stopCluster(cl)
  }else{
    registerDoMC(old)
  }
  # 4. saving iNEXT results
  saveRDS(rare, iNEXT.out)
}else{
  rare <- readRDS(iNEXT.out)
}

# Custom jopning function
my_join <- function(data, by){
  data %>%
    left_join(sample_data(exp) %>% 
                as("data.frame"), by = by)
}

# Buidlign overall table
iNextEst <- map_df(rare, function(d){
  bind_rows(d$iNextEst, .id = "sample.id")
}, .id = "data") %>%
  separateVar() %>%
  my_join(by = "sample.id")

# Building data info table
DataInfo <- map_df(rare, "DataInfo") %>%
  my_join(by = c("site" = "sample.id"))

# Scale x label function
x.lab <- function(x){
  paste0(round(x/1000), "k")
}

# Mean plot (difficult)
# Observed
obs <- iNextEst %>%
  filter(method == "observed") %>%
  group_by(amplicon, env, type.y) %>%
  summarise(lowest = min(m),
            highest = max(m))

# Common sizes (33 samples for biotic, and 15 samples for abiotic)
ms <- iNextEst %>%
  filter(method != "observed") %>%
  group_by(amplicon, type.x, m) %>%
  summarise(c = n()) %>%
  filter(c == 33 | c == 15) %>%
  pull(m)

# Remove observed and filtering only
# common sizes. The mean value of 
# the curve and the confidence interval
# is calculated for each sample type.
#
# Two plots are computed, one for sample
# type and the other for site
c.all <- c("env", "site") %>%
  set_names() %>%
  map(function(i){
    # I have to group dataset by hand
    # since group_by function doesn't
    # seem to work...
    iNextEst.mean <- iNextEst %>%
      filter(method != "observed") %>%
      filter(m %in% ms)
    
    obs <- iNextEst %>%
      filter(method == "observed")
    
    # Grouping by hand... not ideal
    if(i == "env"){
      iNextEst.mean <- iNextEst.mean %>%
        group_by(amplicon, env, m, type.y)
      obs <- obs %>%
        group_by(amplicon, env, type.y)
    }else{
      iNextEst.mean <- iNextEst.mean %>%
        group_by(amplicon, site, m)
      obs <- obs %>%
        group_by(amplicon, site)
    }
    
    # Summarising after grouping
    iNextEst.mean <- iNextEst.mean %>%
      summarise_at(c("qD", "qD.LCL", "qD.UCL"), mean)
    obs <- obs %>%
      summarise(lowest = min(m),
                highest = max(m))
    
    # Final plot
    plots <- iNextEst.mean %>%
      split(.$amplicon) %>%
      map(function(d){
        ggplot(d, aes_string(x = "m", y = "qD", ymin = "qD.LCL", ymax = "qD.UCL",
                             color = i, fill = i, group = i)) +
          ### TO ADD VERTICAL LINES FOR OBSERVED DIVERSITY UNCOMMENT THESE LINES ###
          # geom_vline(data = obs, aes(xintercept = lowest, color = env),
          #            linetype = 2, show.legend = F) +
          # geom_vline(data = obs, aes(xintercept = highest, color = env),
          #            linetype = 3, show.legend = F) +
          geom_ribbon(alpha = .3, color = NA, show.legend = F) +
          geom_line() +
          facet_grid(. ~ type.y, scales = "free",
                     labeller = labeller(type.y = function(x){
                       case_when(
                         x == "biotic" ~ "Crab's organs",
                         T ~ "Environmental samples"
                       )
                     })) +
          scale_color_npg() +
          scale_fill_npg() +
          theme_publication() +
          myTheme +
          xlab("Sequencing depth") +
          ylab("Inverse Simpson index") +
          scale_x_continuous(labels = x.lab,
                             breaks = scales::pretty_breaks(n = 3)) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
          guides(color = guide_legend(nrow = 1, override.aes = list(size = 1)))
      })
    
    if(i == "site"){
      for(n in seq_along(plots)){
        plots[[n]]$facet <- facet_null()
      }
    }
    plots
  })

#### Wilcox test on diversity #####
div <- getDiversityIndexes(otu_table(exp),
                           "col",
                           sample_data(exp)) %>%
  rename_with(~str_remove(.x, "metadata."))

# Formatting table S3
table.s3 <- div %>%
  arrange(env, site) %>%
  select(ID = sample.id, 
         `Sample type` = env, 
         Site = site,
         `Number of ASVs`=Number_of_OTUs,
         Number_of_singletons, Number_of_doubletons,
         `Inverse Simpson index`=InvSimpson,
         `Good's coverage estimator`=Good_coverage_estimator) %>%
  dplyr::rename_all(~str_replace_all(., "_", " "))

goods_range <- table.s3$`Good's coverage estimator` %>%
  range()

# Table 1 main text
# Goods coverage and inverse Simpson 
# delta between observed and extrapolated
table.1 <- c("env", "site") %>% # Once for sample type and once for sites
  set_names() %>%
  map(function(f){
    # Computing mean +/- standard error based on
    # env and sites
    goods <- div %>%
      group_by_at(f) %>%
      summarise(Good_coverage_estimator = sprintf("%.3f ± %.3f", 
                                                  mean(Good_coverage_estimator), 
                                                  sd(Good_coverage_estimator)/sqrt(n())))
    
    # Computing delta between observed and highest extrapolated values
    # for each sample. Values are then grouped by sample type and site
    # and the mean +/- standard error is calculated.
    iNextEst %>%
      group_by_at(c("method", f, "sample.id")) %>%
      summarise(qD = max(qD)) %>%
      filter(method %in% c("extrapolated", "observed")) %>%
      spread(method, qD) %>%
      group_by_at(f) %>%
      mutate(delta = extrapolated - observed) %>%
      summarise(`Inverse Simpson (ext - obs)` = sprintf("%.3f ± %.3f", mean(delta), sd(delta)/sqrt(n()))) %>%
      left_join(goods, by = f)
  }) %>%
  # Final binding and changing names
  bind_rows(.id = "Effect") %>%
  mutate(env = ifelse(is.na(env), 
                      as.character(site), 
                      as.character(env)),
         Effect = case_when(
           Effect == lag(Effect) ~ "",
           T ~ Effect
         )) %>%
  select(-site) %>%
  mutate(Effect = case_when(
    Effect == "env" ~ "Sample type",
    Effect == "site" ~ "Site",
    T ~ Effect
  )) %>%
  dplyr::rename(
    `Good's coverage estimator` = "Good_coverage_estimator",
    `Code` = "env"
  )


# Observed Inverse simpson index for each sample
obs <- iNextEst %>%
  filter(method == "observed")

# Plotting differences based on wilcoxon
# test and corrected using BH p.value correction.
# A plot is drawn for sample types and sites.
all.alpha.plots <- c("env", "site") %>%
  set_names() %>%
  map(function(f){
    # Split contrasts by amplicon
    contrasts <- obs %>%
      split(.$amplicon) %>%
      map(function(x){
        x %>%
          mutate_at(f, ~fct_reorder(., qD, "mean", .desc = T)) %>%
          arrange_at(f)
      })
    
    fct <- get(f, obs)
    col <- pal_npg()(10)[1:nlevels(fct)]
    names(col) <- levels(fct)
    
    # Pairwise wilcoxon test based on sample type
    # or site
    contrasts <- if(f == "env"){
      contrasts %>% 
        map(~pairwise.wilcox.test(.$qD, .$env, p.adjust.method = "BH"))
    }else{
      contrasts %>% 
        map(~pairwise.wilcox.test(.$qD, .$site, p.adjust.method = "BH"))
    }
    
    # Get significant letters
    letters <- contrasts %>%
      map(broom::tidy) %>%
      map(unite, "name", group1, group2, sep = "-") %>%
      map(deframe) %>%
      map(multcompView::multcompLetters) %>%
      map("Letters") %>%
      map_df(enframe, name = f, value = "label", 
             .id = "amplicon")
    
    letters <- obs %>%
      group_by_at(c("amplicon", f)) %>%
      summarise_at(c("qD", "qD.LCL", "qD.UCL"), mean) %>%
      ungroup() %>%
      left_join(letters, by = c("amplicon", f))
    
    # Plot contrasts divided by amplicon type
    c("16S", "ITS") %>%
      set_names() %>%
      map(function(amp){
        cntr <- letters %>%
          filter(amplicon == amp) %>%
          mutate_at(f, ~fct_reorder(., qD, .fun = "mean", .desc = T))
        
        obs %>%
          filter(amplicon == amp) %>%
          mutate_at(f, ~fct_reorder(., qD, .fun = "mean", .desc = T)) %>%
          ggplot(aes_string(x = f, y = "qD", fill = f)) +
          geom_errorbar(data = cntr, aes(y = qD, ymin = qD.LCL, ymax = qD.UCL),
                        width = .3, size = .25, show.legend = F) +
          stat_summary(geom = "col", fun = mean, show.legend = F) +
          geom_text(data = cntr, aes(y = qD.UCL, label = label), 
                    vjust = -1, size = myTheme$text$size / .pt) +
          myTheme +
          scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
          scale_fill_manual(values = col) +
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank())
      })
  })

# Saving the final plot (Figure 2):
#  panels a and c : difference based on 16S counts
#  panels b and d : difference based on ITS counts
layout <- "AAAABB
           CCCCDD"
figure.2 <- c.all$env$`16S` + all.alpha.plots$env$`16S` + 
  c.all$env$`ITS` + all.alpha.plots$env$`ITS` +
  plot_layout(guides = "collect", design = layout) +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.tag = element_text(vjust = -3, face = "bold", 
                                size = myTheme$text$size + 1),
        plot.margin = unit(c(0,.1,0,.1), "lines"),
        plot.background = element_blank())

# Figure S3
figure.s3 <- c.all$site$`16S` + all.alpha.plots$site$`16S` + 
  c.all$site$`ITS` + all.alpha.plots$site$`ITS` +
  plot_layout(guides = "collect", design = layout) +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.tag = element_text(vjust = -3, face = "bold", 
                                size = myTheme$text$size + 1),
        plot.margin = unit(c(0,.1,0,.1), "lines"),
        plot.background = element_blank())


# ---- betaDiversity ------------------------------------------

# Splitting experiment according to "superdom" (16S and ITS) and
# environment type (biotic or abiotic). Each dataset is then
# normalized using relative abundance
exp.split <- splitBy(exp, taxa = "superdom", sample = "type",
                     both = F, keep.zeroes = F)

# If beta-diversity is quantitative remove 
# datatset with both 16S and ITS
if(!binary){
  extr <- grep("16S\\|ITS", names(exp.split))
  if(length(extr) > 0)
    exp.split <- exp.split[-extr]
}

# Transform into relative abundace
exp.rel <- exp.split %>% 
  map(transform_sample_counts, function(x) (x / sum(x)))

# permanova using crab as blocking factor
# for biotic factors only
permanova <- exp.rel %>%
  map(function(e){
    x <- otu_table(e) %>%
      as("matrix")
    
    d <- sample_data(e) %>%
      as("data.frame") %>%
      mutate(crab.mod = paste0(crab, "_", site))
    
    design <- if(all(d$type == "biotic")){
      x ~ env + site + crab.mod
    }else{
      x ~ env + site
    }
    set.seed(123)
    adonis2(design, d, 1000, "bray", binary = binary)
  }) %>%
  map_df(as_tibble, rownames = "Effect", 
         .id = "data") %>%
  separateVar()

# Betadisper (calculated on biotic and abiotic factors)
b.disp <- exp.rel %>%
  map(function(e){
    x <- otu_table(e) %>%
      as("matrix")
    
    grp <- c("env", "site")
    
    set.seed(123)
    grp %>%
      set_names() %>%
      map(get_variable, physeq = e) %>%
      map(~betadisper(d = vegdist(x, "bray", binary = binary), 
                      group = .))
  }) 

# ANOVA table for beta dispersion
table.s5 <- b.disp %>%
  map_df(function(x) {
    map(x, anova) %>%
      map_df(broom::tidy, .id = "Effect")
  }, .id = "data") %>%
  separateVar() %>%
  filter(term == "Groups") %>%
  select(-term) %>%
  mutate(type = case_when(
    type == "biotic" ~ "Crab's organs",
    type == "abiotic" ~ "Environmental samples",
  )) %>%
  mutate(Effect = case_when(
    Effect == "env" ~ "Sample type",
    Effect == "site" ~ "Site",
    T ~ Effect
  )) %>%
  mutate(p.value = sprintf("%.3f%s", p.value, 
                           symnum(p.value, c(0,0.01,0.05,1), c("**", "*", "")))) %>%
  select(amplicon, type, Effect, p.value) %>%
  spread(Effect, p.value) %>%
  mutate(amplicon = case_when(
    amplicon == lag(amplicon) ~ "",
    T ~ amplicon)) %>%
  rename_with(.fn = str_to_title)

##### PCoA on splitted data #####
# PCoA using Bray-Curtis distance
set.seed(1783)
ords <- exp.rel %>%
  map(otu_table) %>%
  map(as, "matrix") %>%
  map(vegdist, method = "bray") %>%
  map(cmdscale, eig = T)

# Get ordination data using the custom 
# "getOrdination" funciton
ords.data <- exp.rel %>%
  map2(ords, function(e, o) {
    meta <- sample_data(e) %>% as("data.frame")
    getOrdination(o, meta = meta, centroids = c("env", "site"), 
                  spider = T, hull = T, sample = "sample.id")
  })


# Function for plotting data according to
# a given group.
plot_pcoa <- function(pcoas, grp){
  main.data <- ords.data %>%
    map_df("data", .id = "data") %>%
    separateVar()
  
  spiders <- ords.data %>%
    map("spider") %>%
    map_df(grp, .id = "data") %>%
    separateVar()
  
  centroids <- ords.data %>%
    map("centroids") %>%
    map_df(grp, .id = "data") %>%
    separateVar()
  
  variance <- ords.data %>%
    map_df("variance", .id = "data") %>%
    separateVar()

  
  variance <- variance %>%
    summarise(x = range(x),
              y = range(y))
  labs <- sprintf("PCoA%d (%.2f - %.2f%%)", 1:2, variance[1,], variance[2,])
  ggplot(main.data, aes_string(x = "x", y = "y", color = grp)) +
    geom_hline(yintercept = 0, linetype = 3, size = .25) +
    geom_vline(xintercept = 0, linetype = 3, size = .25) +
    geom_segment(data = spiders, aes(xend = xend, yend = yend),
                 size = .2, linetype = 1, show.legend = F) +
    geom_point(size = 1) +
    # geom_label(data=centroids, aes_string(label = grp), size = 2,
    #            label.padding = unit(.1, "lines"), show.legend = F) +
    geom_point(data=centroids, shape = 21, fill = "white", show.legend = F,
               size = 1.5, stroke = .2) +
    facet_grid(amplicon ~ type, scales = "free",
               labeller = labeller(amplicon = labelFrmt,
                                   type = labelFrmt)) +
    scale_color_npg() +
    theme_bw(base_size = 6, base_family = "sans") +
    theme_publication(grid = F) +
    xlab(labs[1]) +
    ylab(labs[2]) +
    guides(color = guide_legend(nrow = 1,
                                override.aes = list(size = 3))) +
    theme(legend.title = element_blank(),
          panel.spacing.x = unit(1, "lines"),
          panel.spacing.y = unit(.7, "lines"))
}

# PCoA main text
a <- plot_pcoa(ords.data, "env")
# PCoA SM
a.supp <- plot_pcoa(ords.data, "site")


# Better to make a table with results instead of a plot
table.s4 <- permanova %>%
  mutate(Effect = case_when(
    Effect == "env" ~ "Sample type",
    Effect == "site" ~ "Site",
    Effect == "crab.mod" ~ "Crab",
    T ~ Effect
  )) %>%
  mutate(type = case_when(
    type == "biotic" ~ "Crab's organs",
    type == "abiotic" ~ "Environmental samples",
    T ~ "Both"
  )) %>%
  mutate(sign = symnum(`Pr(>F)`, c(0, 0.01, 0.05, 1), 
                       c("**", "*", ""), na = F)) %>%
  filter(Effect != "Total") %>%
  mutate(R2 = sprintf("%.3f%s", R2, sign)) %>%
  select(amplicon, type, Effect, R2) %>%
  spread(Effect, R2, fill = "-") %>%
  mutate(amplicon = case_when(
    amplicon == lag(amplicon) ~ "",
    T ~ amplicon
  )) %>%
  rename_with(str_to_title) %>%
  relocate(Residual, .after = last_col())

##### Pairwise PERMANOVA #####
# Pairwise permanova on env and site
permanova.pair <- exp.rel %>%
  map_df(function(e){
    x <- otu_table(e) %>%
      as("matrix")
    
    d <- sample_data(e) %>%
      as("data.frame") %>%
      mutate(crab.mod = paste0(crab, "_", site))
    
    design <- if(all(d$type == "biotic")){
      x ~ env + site + crab.mod
    }else{
      x ~ env + site
    }
    set.seed(123)
    c("env", "site") %>%
      set_names() %>%
      map_df(~pairwise.adonis(x, pull(d, .), "bray", "BH", binary = binary),
             .id = "Effect")
    
  }, .id = "data") %>%
  separateVar() %>%
  separate(pairs, into = c("group1", "group2"), sep = " vs ")

# Pairwise R2 plot
b <- permanova.pair %>%
  mutate(R2 = ifelse(p.adjusted < 0.05, R2, NA),
         group1.b = ifelse(amplicon == "ITS", group2, group1),
         group2 = ifelse(amplicon == "ITS", group1, group2)) %>%
  select(-group1) %>%
  dplyr::rename(group1=group1.b) %>%
  split(.$Effect) %>%
  map(function(d){
    max <- ifelse("GI" %in% pull(d, group1), 4, 3)
    ggplot(d, aes(x = group1, y = group2, fill = R2)) +
      geom_point(shape = 21, aes(size = R2)) +
      annotate(geom = "segment", x = 1, y = 1, xend = max, yend = max,
               lineend = "round") +
      annotate(x = max - .3, y = max, geom = "text", label = "16S", 
               size = 2, hjust = .5, fontface='bold', angle = 47) +
      annotate(x = 1.3, y = 1, geom = "text", label = "ITS", 
               size = 2, hjust = .5, fontface='bold', angle = 47) +
      facet_wrap(~ type, scales = "free", labeller = labeller(type = labelFrmt)) +
      myTheme +
      coord_cartesian(clip = "off") +
      scale_fill_distiller(palette = "RdYlBu") +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "right",
            legend.direction = "vertical",
            aspect.ratio = 1) +
      guides(fill=guide_legend(title = bquote(PERMANOVA~R^2)), 
             size = guide_legend(title = bquote(PERMANOVA~R^2)))
  })

# Picking plots
b.supp <- b$site
b <- b$env

# Saving using cowplot to avoid unwanted alignements
# generate dby patchwork
lab.size <- myTheme$text$size + 1
figure.3 <- plot_grid(a, b, ncol = 1, rel_heights = c(2,1), 
                      labels = c("a", "b"), label_size = lab.size)
figure.s4 <- plot_grid(a.supp, b.supp, ncol = 1, rel_heights = c(2,1), 
                           labels = c("a", "b"), label_size = lab.size)

# ---- commonASVs ------------------------------

###### Sharing across groups #######
# Buidling datasets (one for 16S an one for ITS)
exp.upset <- splitBy(exp, taxa = "superdom", 
        both = F, keep.zeroes = F)

# Making upset plots
upsets <- exp.upset %>%
  map(function(e){
    # Intersections of 2nd degree only (for 
    # sample type)
    inters <- sample_data(e) %>%
      pull(env) %>% unique() %>%
      as.character() %>%
      (function(x) {
        res <- combn(x, m = 2, simplify = F)
        append(res, x)
      })
    
    # Building dataset for sample type
    byEnv <- merge_samples(e, "env") %>%
      otu_table() %>% t() %>% sign() %>%
      as("matrix") %>% data.frame()
    
    # Building dataset for sites
    bySites <- merge_samples(e, "site") %>%
      otu_table() %>% t() %>% sign() %>%
      as("matrix") %>% data.frame()
    
    list(
      bySites = upset(bySites, nsets = 3, show.numbers = "no"),
      byEnv = upset(byEnv, nsets = 8, intersections = inters,
                    show.numbers = "no")
    )
  })

# Building final plot and saving
figure.s5 <- upsets %>%
  map2(seq_along(upsets), function(x, i){
    plots <- map(x, function(y){
      plot_grid(y$Main_bar, y$Matrix, ncol = 1,
                align = "hv", axis = "b")
    })
    labs <- if(i == 1){
      c("a", "b")
    }else{
      c("c", "d")}
    plot_grid(plotlist = plots, ncol = 2,
              rel_widths = c(1, 3),
              label_size = myTheme$text$size + 1,
              labels = labs)
  }) %>%
  plot_grid(plotlist = ., ncol = 1)