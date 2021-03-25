# ---- formatData --------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(scales)
  library(ggbeeswarm)
  library(Biostrings)
})
source("scripts/utils.R")

## 33 gills  (11*3 sites)
## 33 indgut (11*3 sites)
## 33 midgut (11*3 sites)
## 33 gonads (11*3 sites)
## 15 soils  (5*3 sites)
## 15 litter (5*3 sites)
## 15 water  (5*3 sites)
## 15 wd     (5*3 sites)
## tot:
totSamples <- (33 * 4) + (15 * 4)

## CH = GI (gills)
## HG = HIN (indgut)
## MG = MID (midgut)
## (XX)L = litter
## (XX)S = soil
## (XX)W = water
## (XX)WD = water debris

# Counts and taxa
counts <- list.files(pattern = "all_seqtab_nochim.rds", 
                     recursive = T, full.names = T)
taxa <- list.files(pattern = "tax_assignments.rds",
                   recursive = T, full.names = T)


# Trackback
track <- list.files(pattern = "track_back.tsv",
                    recursive = T, full.names = T)
track <- map(track, read.table, sep = "\t", 
             header = T, row.names = 1) %>%
  set_names(str_extract(counts, "16s|its"))

getMetadata <- function(id){
  # Metadata are included in the file IDs following
  # a shared scheme with some exceptions:
  #
  # 1. GI (gills) is some samples has been reported as CH
  # 2. HG (hidgut) is some samples has been reported as HIN
  # 3. MG (imdgut) is some samples has been reported as MID
  # 4. ITS samples have a '-' before sample number
  # 5. Some replicates have been flagged with a '-'
  #    followed by a number
  # 6. Environmental samples have the samplig site reported
  #    in their ID whereas crab samples are reported in
  #    chunck of 11 samples for each site
  
  # Step 1: change CH, HIN e MID to GI, HG, MG.
  #         Remove '-' before sample number
  #         Remove '-2' from replicates
  idStn <- id %>%
    str_replace("CH", "GI") %>%
    str_replace("HIN", "HG") %>%
    str_replace("MID", "MG") %>%
    str_remove("(?<=[0-9])-[0-9]") %>%
    str_replace("-", "") %>%
    str_extract(".*(?=_)")
  
  # Step 2: get label and sample number
  lab <- str_extract(idStn, "[A-Z]+(?=[0-9])")
  ord <- as.numeric(str_extract(idStn, "(?<=[A-Z])[0-9]+"))
  
  # Step 3: split environmental samples from crab samples
  #         (according to the description provided in 6)
  envs <- "(L|S|W|WD)$"
  env <- str_extract(lab, envs)
  
  # Step 4: Set "unknown" site for crab samples
  type <- ifelse(is.na(env), "biotic", "abiotic")
  site <- ifelse(is.na(env), "unkn", str_replace(lab, envs, ""))
  env <- ifelse(is.na(env), lab, env)
  
  # Step 6: add site to crab samples
  sites <- c("TKP", "SK", "SH")
  
  res <- data.frame(id, idStn, lab, ord, env, type, site,
                    stringsAsFactors = F) %>%
    # grouping chuncks
    arrange(env, site, ord) %>%
    group_by(env, site) %>%
    # ranking
    mutate(ord = dense_rank(ord)) %>%
    # find group
    mutate(n = ceiling(ord/11)) %>%
    ungroup() %>%
    # add label
    mutate(site = ifelse(site == "unkn", sites[n], site)) %>%
    # final format
    group_by(env, site) %>%
    mutate(sample = dense_rank(ord),
           sample.id = paste(env,site,sample, sep = "_"),
           crab = ifelse(type == "biotic", sample, 0)) %>%
    select(id, idStn, sample.id, env, site, sample, crab, type)

  # Reordering
  ord <- res %>% pull(id) %>%
    match(id, .)
  
  res[ord,]
}

# Building final datatsets
f <- function(d){
  m <- getMetadata(rownames(d))
  data.frame(m, d, row.names = NULL)
}
track <- track %>%
  map_df(f, .id = "gene")

# Reading tables (counts + taxonomies)
x <- map(counts, readRDS) %>%
  set_names(str_extract(counts, "16s|its"))

tax <-  map(taxa, readRDS) %>%
  set_names(str_extract(taxa, "16s|its"))

# Building raw phyloseq objects
phylos <- map2(x, tax, function(a, b){
  # metadata from rownames
  m <- getMetadata(rownames(a))
  # order of samples
  ord <- with(m, order(env, site, sample))
  
  # reordering
  m <- m[ord,]
  a <- a[ord,]
  
  # naming
  # n <- paste0("ASV_", 1:ncol(a))
  # colnames(a) <- n
  rownames(a) <- NULL
  
  # rownames(b) <- n
  
  # phyloseq object
  phyloseq(otu_table(a, taxa_are_rows = F),
           sample_data(m),
           tax_table(b))
})
# Removing T sample from 16S
phylos$`16s` <- subset_samples(phylos$`16s`, env != "T")

# Removing contaminants
getContaminants <- function(p){
  # known taxa (domain different from NA)
  known <- subset_taxa(p, !is.na(domain))
  # removing choloroplast and mitochondria
  no_chloroplast <- subset_taxa(known, order != "Chloroplast")
  no_mitochondria <- subset_taxa(no_chloroplast, order != "Mitochondria")
  
  # Tabling data
  m <- sample_data(p) %>%
    as("data.frame")
  
  res <- data.frame(m,
             known = sample_sums(known),
             no_chloroplast = sample_sums(no_chloroplast),
             no_mitochondria = sample_sums(no_mitochondria))
  list(data = res, phylo = no_mitochondria)
}

# Apply function
res <- map(phylos, getContaminants)

# Get phyloseq objects without contaminants
phylos <- map(res, "phylo")

# Getting number of reads and number of
# asv per sample
data <- map_df(phylos, function(p){
  n.asv <- otu_table(p) %>%
    as("matrix") %>%
    apply(1, function(x) sum(x > 0))
  n.reads <- sample_sums(p)
  d <- sample_data(p) %>%
    as("data.frame")
  cbind(d, n.asv, n.reads)
}, .id = "gene")

# Testing if the number of sample is correct
ok <- data %>%
  group_by(gene, env, sample.id) %>%
  tally() %>%
  tally() %>%
  pull(n) %>%
  sum() %>%
  "=="(totSamples * 2)
if(!ok){
  warning("The number of samples is not correct")
}

# ---- correlations --------------------------------
# Tecnical replicates:
reps.sample <- c("GI_SH_4", "GI_SH_9", "GO_SH_8", "GO_TKP_7", 
                 "HG_SK_6", "HG_SH_1", "MG_SK_11", "MG_SH_3")

# Saving table ordered by sample.id
data %>%
  arrange(gene, sample.id) %>%
  mutate(replicate = case_when(
           sample.id %in% reps.sample & 
             gene == "its" ~ "R",
           TRUE ~ "")) %>%
  group_by(sample.id, gene) %>%
  mutate(replicate = case_when(
    replicate == "R" ~ paste0(replicate, 1:n()),
    TRUE ~ "")) %>%
  write_excel_csv2("reads_asv.csv")

# Selecting samples
reps <- phylos$its %>%
  subset_samples(sample.id %in% reps.sample) %>%
  transform_sample_counts(function(x) (x / sum(x)) * 1000) %>%
  filter_taxa(function(x) sum(x) > 0, T)

m <- otu_table(reps) %>%
  as("matrix") %>%
  t()
s <- sample_data(reps) %>%
  as("data.frame") %>%
  mutate(sample.name = sample_names(reps))

# Spearman's correlation
correlation <- data.frame(m, check.names = F, row.names = NULL) %>%
  pairwise(all = F, diag = F, join.var = s, 
           join.by = "sample.name", fun = function(x) log2(x + 1)) %>%
  filter(sample.id.x == sample.id.y)


correlate <- function(x, y, method = "spearman"){
  extr <- x > 0 & y > 0
  cor(x[extr], y[extr], method = method)
}

# Accuracy
accuracy <- function(x, y){
  tp <- sum(x > 0 & y > 0)
  tn <- sum(x == 0 & y == 0)
  
  fpfn <- sum( (x > 0 & y == 0) | (x == 0 & y > 0) )
  (tp + tn)/(tp + tn + fpfn)
}

# Get Spearman's rho
rhos <- correlation %>%
  group_by(sample.id.x, env.x) %>%
  summarise(N = n(),
            r = correlate(x, y),
            tp = sum(x > 0 & y > 0),
            tn = sum(x == 0 & y == 0),
            fpfn = sum( (x > 0 & y == 0) | (x == 0 & y > 0) ),
            accuracy = (tp + tn)/(tp + tn + fpfn)) %>%
  setNames(c("Sample id", "env", "N", "rho", "TP", "TN", "FP + FN", "Accuracy")) %>%
  ungroup()

# Plotting
labs <- rhos %>%
  group_by(env) %>%
  summarize(lab = sprintf("%.2f", mean(rho))) %>%
  mutate(x = 1,
         y = max(correlation$y)) %>%
  dplyr::rename(env.x = env) 

# Figure S2
figure.s2 <- correlation %>%
  sample_frac() %>%
  filter(x > 0, y > 0) %>%
  ggplot(aes(x = x, y = y, color = sample.id.x)) +
  geom_point(show.legend = F, size = pt2ggSize(3)) +
  # Adding rho means
  geom_text(data = labs, aes(label = lab), color = "black",
            show.legend = F, hjust=0, vjust = 1,
            size = pt2ggSize(6)) +
  facet_wrap(~ env.x) +
  # Coord equals
  coord_equal(xlim = c(0, 10), ylim = c(0,10)) +
  # Ab lines for perfect correlation
  geom_abline(slope = 1, intercept = 0, linetype = 3,
              size = pt2ggSize(1)) +
  # Graphical parameters
  theme_bw(base_size = 8, base_family = "sans") +
  theme_publication(grid = F) +
  scale_color_brewer(palette = "Paired") +
  xlab("Abundace of SVs in replicate 1") +
  ylab("Abundace of SVs in replicate 2")


# ---- trackback -----------------------------
data <- map_df(res, "data", .id = "gene")
track <- track %>% 
  right_join(data) %>%
  mutate(sample = as.character(sample),
         crab = as.character(crab)) %>%
  select(-id) %>%
  group_by_if(is.character) %>%
  summarise_if(is.numeric, sum) %>%
  dplyr::select(-derepR, -denoisedR) %>%
  dplyr::rename(`Raw reads` = "raw",
                `PCR primers removal` = "no.adapt",
                `Quality filtering` = "filtered",
                `Dereplication` = "derepF",
                `Denoising` = "denoisedF",
                `Merging` = "merged",
                `Chimera removal` = "no.chim",
                `Known Domain` = "known",
                `Chloroplast removal` = "no_chloroplast",
                `Mitochondria removal` = "no_mitochondria")


track.back <- track %>%
  ungroup() %>%
  gather("step", "count", 
         -type, -sample, -crab, -env, -site, -sample.id, -idStn, -gene) %>%
  mutate(step = fct_inorder(step),
         env = fct_relevel(env, "HG", "MG", "GI", "GO")) 

track.back <- track.back %>%
  filter(!(step %in% c("Chloroplast removal", "Mitochondria removal") & gene == "its")) %>%
  mutate(gene = case_when(
    gene == "16s" ~ "16S rRNA gene (V3 - V4)",
    T ~ "Internal transcribed spacer (ITS1)"
  )) %>%
  group_by(type) %>%
  mutate(type = case_when(
    type == "biotic" ~ paste0("Crab tissues\n(", length(unique(sample.id)), " samples)"),
    T ~ paste0("Envrironmental samples\n(", length(unique(sample.id)), " samples)")
  ))

its.removed <- track.back %>%
  filter(count == 0)

figure.s1 <- track.back %>%
  filter(count > 0) %>%
  ggplot(aes(x = step, y = count + 1)) +
  geom_boxplot(outlier.shape = NA, size = .25) +
  geom_point(data = its.removed, aes(group = sample.id),
             alpha = .4, size = .8, 
             position = position_dodge2(width = .4)) +
  geom_text(data = its.removed, aes(group = sample.id, label = sample.id),
             size = 4 / .pt, 
             position = position_dodge2(width = 1.6),
            hjust = .2, angle = 90) +
  # geom_quasirandom(aes(color = type), alpha = .4, size = .5) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                labels = function(x) {
                  ifelse(x == 1, 0, 
                         trans_format("log10", math_format(10^.x))(x))
                }) +
  facet_grid(type ~ gene, space = "free_x", scales = "free_x") +
  theme_bw(base_family = "sans", base_size = 8) +
  theme_publication(grid = "y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  ylab("Number of reads retained (log)") +
  xlab("Amplicon variant detection step") +
  coord_cartesian(ylim = c(1, 1e6))



# Percentage of retained reads in respect to
# previous step
perc.retained <- track.back %>%
  group_by(gene, step) %>%
  summarise(retained = sum(count)) %>%
  mutate(retained = (retained / lag(retained)) * 100) %>%
  ungroup() %>%
  split(.$gene) %>%
  map(select, -gene) %>%
  map(deframe)

# ---- saveData ---------------------------------------
# Merging data based on sample.id
phylos <- phylos %>%
  map(merge_samples_mod, group = "sample.id")

# saving sequences (tree eventually) and
# change ASV name
phylos <- phylos %>%
  imap(function(p, i){
    seq <- DNAStringSet(taxa_names(p))
    n <- paste0("ASV_", i, "_", 1:ntaxa(p))
    names(seq) <- n
    taxa_names(p) <- n
    
    out <- paste0("output/", i, "_asv.fasta")
    writeXStringSet(seq, out)
    
    p
  })


# Saving merged experiment
exp <- merge_phyloseq(phylos$`16s`, phylos$its)

# Refactoring
sample_data(exp) <- sample_data(exp) %>%
  as("data.frame") %>%
  #  rownames_to_column("sample.id") %>%
  changeLevels() %>%
  data.frame(row.names = .$sample.id) %>%
  sample_data()

taxa <- tax_table(exp) %>% as("matrix")
superdom <- ifelse(taxa[,"domain"] == "Fungi", "ITS", "16S")
tax_table(exp) <- tax_table(cbind(superdom, taxa))

saveRDS(exp, "./output/exp.rds")