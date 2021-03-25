# Converting ggplot2 point size to standard size
pt2ggSize <- function(pt) {
  # The number of points corresponding to size 1 in ggplot
  magic_number <- 2.14
  pt/magic_number
}

# Overall theme set function
theme_publication <- function(lines = .5, rects = .5, text = 8, ticks.len = 2, legend = 10,
                              legend.position = "bottom", legend.direction = "horizontal",
                              grid = T, ticks = T, axis.line = F) {
  tm <- theme(text = element_text(size = text),
              plot.title = element_text(face = "bold",
                                        size = rel(1.2), hjust = 0.5),
              line = element_line(size = pt2ggSize(lines), lineend = "square"),
              rect = element_rect(size = pt2ggSize(rects)),
              
              # panel.background = element_rect(colour = NA),
              # plot.background = element_rect(colour = NA),
              # panel.border = element_rect(colour = NA),
              
              axis.title = element_text(size = rel(1)),
              axis.title.y = element_text(angle=90, vjust = 2, size = rel(1)),
              axis.title.x = element_text(vjust = -0.2, size = rel(1)),
              
              axis.text = element_text(size = rel(1)), 
              
              # axis.line = element_line(colour="black"),
              axis.ticks = if(ticks) element_line() else element_blank(),
              axis.ticks.length = unit(ticks.len, "pt"),
              
              panel.grid.major.x = if(grid == "x" | grid == T) element_line(colour="#f0f0f0") else element_blank(),
              panel.grid.major.y = if(grid == "y" | grid == T) element_line(colour="#f0f0f0") else element_blank(),
              panel.grid.minor = element_blank(),
              
              legend.key = element_rect(colour = NA),
              legend.position = legend.position,
              legend.direction = legend.direction,
              legend.key.size= unit(legend, "pt"),
              legend.margin = margin(0, 0, 0, 0, "cm"),
              
              legend.title = element_text(face="bold", size = rel(.8)),
              legend.text = element_text(size = rel(.8)),
              
              plot.margin=unit(c(3,2,3,3),"mm"),
              
              strip.background=element_blank(),
              strip.text = element_text(size = rel(1))
  )
  
  if(axis.line){
    tm + theme(axis.line = element_line(size = pt2ggSize(lines), lineend = "square"))
  }else{
    tm
  }
}


# Defining theme
myTheme <- theme_bw(base_size = 8, base_family = "sans") +
  theme_publication() +
  theme(strip.background = element_blank())

# Changing levels according to a given order
changeLevels <- function(data){
  data %>%
    mutate(crab = factor(crab),
           env = fct_relevel(env, "MG", "HG", "GO", "GI",
                             "L", "S", "W", "WD"),
           site = factor(site),
           type = ifelse(env %in% c("L", "S", "W", "WD"),
                         "abiotic", "biotic"))
}

# GGplot theme
theme_publication <- function(lines = .5, rects = .5, text = 8, ticks.len = 2, legend = 10,
                              legend.position = "bottom", legend.direction = "horizontal",
                              grid = T, ticks = T, axis.line = F) {
  tm <- theme(text = element_text(size = text),
              plot.title = element_text(face = "bold",
                                        size = rel(1.2), hjust = 0.5),
              line = element_line(size = pt2ggSize(lines), lineend = "square"),
              rect = element_rect(size = pt2ggSize(rects)),
              
              # panel.background = element_rect(colour = NA),
              # plot.background = element_rect(colour = NA),
              # panel.border = element_rect(colour = NA),
              
              axis.title = element_text(size = rel(1)),
              axis.title.y = element_text(angle=90, vjust = 2, size = rel(1)),
              axis.title.x = element_text(vjust = -0.2, size = rel(1)),
              
              axis.text = element_text(size = rel(1)), 
              
              # axis.line = element_line(colour="black"),
              axis.ticks = if(ticks) element_line() else element_blank(),
              axis.ticks.length = unit(ticks.len, "pt"),
              
              panel.grid.major.x = if(grid == "x" | grid == T) element_line(colour="#f0f0f0") else element_blank(),
              panel.grid.major.y = if(grid == "y" | grid == T) element_line(colour="#f0f0f0") else element_blank(),
              panel.grid.minor = element_blank(),
              
              legend.key = element_rect(colour = NA),
              legend.position = legend.position,
              legend.direction = legend.direction,
              legend.key.size= unit(legend, "pt"),
              legend.margin = margin(0, 0, 0, 0, "cm"),
              
              legend.title = element_text(face="bold", size = rel(.8)),
              legend.text = element_text(size = rel(.8)),
              
              plot.margin=unit(c(3,2,3,3),"mm"),
              
              strip.background=element_blank(),
              strip.text = element_text(size = rel(1))
  )
  
  if(axis.line){
    tm + theme(axis.line = element_line(size = pt2ggSize(lines), lineend = "square"))
  }else{
    tm
  }
}

pt2ggSize <- function(pt) {
  # The number of points corresponding to size 1 in ggplot
  magic_number <- 2.14
  pt/magic_number
}


splitBy <- function(data, sample = NULL, taxa = NULL, both = F, keep.zeroes = F){
  # Splits a phyloseq object into chunks depending on
  # taxa rank and sample variable selected.
  # If both == T then variable are combined using the "|"
  # symbol to select all combinations
  library(tidyverse)
  library(phyloseq)
  
  if(sum(is.null(c(sample, taxa))) == 2){
    return(data)
  }
  
  taxa.fct <- ""
  if(!is.null(taxa)){
    taxa.fct <- as(tax_table(data), "matrix")[,taxa] %>%
      as.vector() %>% unique()
  }
  
  sample.fct <- ""
  if(!is.null(sample)){
    sample.fct <- get_variable(data, sample) %>%
      as.character() %>% unique()
  }
  
  if(both){
    taxa.fct <- append(taxa.fct, paste0(taxa.fct, collapse = "|"))
    sample.fct <- append(sample.fct, paste0(sample.fct, collapse = "|"))
  }
  
  fct <- expand.grid(taxa.fct, sample.fct, stringsAsFactors = F) %>%
    distinct(Var1, Var2)
  
  pruneTaxa <- function(exp, taxa, subset, grep = F){
    t <- tax_table(exp) %>%
      as("matrix") %>%
      as_tibble() %>%
      pull(taxa)
    
    extr <- if(grep){
      grepl(subset, t)
    }else{
      t == subset
    }
    prune_taxa(extr, exp)
  }
  
  pruneSamples <- function(exp, sample, subset, grep = F){
    s <- get_variable(exp, sample)
    
    extr <- if(grep){
      grepl(subset, s)
    }else{
      s == subset
    }
    prune_samples(extr, exp)
  }
  
  l <- as.vector(t(fct[3,]))
  res <- apply(fct, 1, function(l){
    sub <- data
    
    if(!is.null(taxa)){
      g <- grepl("\\|", l[1])
      sub <- pruneTaxa(sub, taxa, l[1], grep = g)
      
      sub.data <- sample_data(sub) %>%
        as("data.frame") %>%
        cbind(taxa.cntr = l[1])
      
      sample_data(sub) <- sample_data(sub.data)
    }
    
    if(!is.null(sample)){
      g <- grepl("\\|", l[2])
      sub <- pruneSamples(sub, sample, l[2], grep = g)
      
      sub.data <- sample_data(sub) %>%
        as("data.frame") %>%
        cbind(sample.cntr = l[2])
      
      sample_data(sub) <- sample_data(sub.data)
    }
    
    if(!keep.zeroes){
      sub <- prune_samples(sample_sums(sub) > 0, sub)
      filter_taxa(sub, function(x) sum(x) > 0, TRUE)
    }else{
      sub
    }
  })
  
  names(res) <- apply(fct, 1, paste0, collapse = "__") %>%
    str_remove("__$|^__")
  res
}


pairwise <- function(data, ..., all = F, diag = F, fun = NULL,
                     join.var = NULL, join.by = NULL){
  #' This function take a data.frame as input and
  #' computes all combinations of variables.
  #' A molteni data.frame will be produced with
  #' at least 5 columns:
  #'
  #' x:          the colum with the first values
  #' y:          the column with the second values
  #' group.x:    the name of the variable in x
  #' group.y:    the name of the variable in y
  #' group.both: the contrast itself
  #'
  #' If specified, the function 'fun' will be
  #' applyed to all variables selected (or
  #' not excluded).
  #'
  #' args:
  #'
  #' data:     the data.frame to transform
  #' all:      if TRUE all contrast will be reported
  #' diag:     if TRUE contrast between the same variable
  #'           will be reported
  #' fun:      the function that will be applied to each
  #'           variable
  #' join.var: if variables are defined in another data frame
  #'           it can be joined to the results produced by
  #'           this function.
  #' join.by:  if specified the join.var data frame will be
  #'           joined by the column selected.
  #' ...:      variables to include/exclude from the contrasts.
  #'           Variables excluded will be added as new column
  #'           and can be used for additional grouping (see
  #'           tidyr::gather function for a similar syntax)
  #'
  #' returns:
  #' A molteni data.frame with al least 5 columns
  
  # Required libraries
  library(tidyr)
  library(dplyr)
  library(forcats)
  library(purrr)
  
  # Parsing args
  quo_args <- quos(...)
  if(length(quo_args) == 0){ # If no variables are
    m <- data                # provided keep the whole
    f <- data.frame()        # data frame
  }else{
    m <- data %>% select(!!!quo_args)
    f <- data %>%
      select(setdiff(colnames(data), colnames(m)))
  }
  
  # Make contrasts
  c <- expand.grid(colnames(m),
                   colnames(m))
  colnames(c) <- c("key1", "key2")
  
  # If fun is provided transform
  # all variable selected using fun
  if(!is.null(fun))
    m <- m %>% mutate_all(fun)
  
  # Filter contrasts based on options
  if(all){ # If all keep all contrast
    if(!diag) # If diag = FALSE remove diags
      c <- c %>% filter(key1 != key2)
  }else{ # If all = FALSE remove redundant contrasts
    d <- c %>% filter(key1 == key2)
    c <- data.frame(t(combn(colnames(m), 2)))
    colnames(c) <- c("key1", "key2")
    if(diag){ # Same as above
      c <- suppressWarnings(bind_rows(c, d))
    }
    rm(d)
  }
  
  # Order levels based on occurrence
  ord <- levels(fct_infreq(c$key1))
  
  # Make the final data
  # c <- c %>% mutate_all(as.character) %>%
  #   rowwise() %>%
  #   do( # repeat all variables as x and y
  #     data.frame(x = m[,.$key1],
  #                y = m[,.$key2],
  #                group.x = .$key1,
  #                group.y = .$key2,
  #                group.both = paste(.$key1, .$key2, sep = "_vs_"),
  #                stringsAsFactors = F) %>%
  #       bind_cols(f)
  #   ) %>% ungroup()
  
  c <- c %>% mutate_all(as.character) %>%
    pmap(function(key1, key2){
      res <- data.frame(x = m[,key1],
                        y = m[,key2],
                        group.x = key1,
                        group.y = key2,
                        group.both = paste(key1, key2, sep = "_vs_"),
                        stringsAsFactors = F)
      if(nrow(f) > 0){
        res <- cbind(res,f)
      }
      return(res)
    }) %>% bind_rows() %>% ungroup()
  
  # Joining variables if specified
  if(!is.null(join.var)){
    join.var <- as.data.frame(join.var)
    if(is.null(join.by)){
      library(tibble)
      join.var <- join.var %>%
        rownames_to_column(var = "joining.by")
      join.by <- "joining.by"
    }
    c <- c %>% left_join(join.var, by = c("group.x" = join.by)) %>%
      left_join(join.var, by = c("group.y" = join.by))
  }
  
  ord.x <- ord[ord %in% unique(c$group.x)]
  ord.y <- ord[ord %in% unique(c$group.y)]
  
  c <- c %>% mutate(group.x = fct_relevel(group.x, ord.x),
                    group.y = fct_relevel(group.y, ord.y))
  
  return(c)
}

#' Merging phyloseq samples by factor
#'
#' Modified  version of merge_samples function
#' from phyloseq package. This function will
#' handle sample data correctly.
#'
#' @param phylo a phyloseq object
#' @param group the name of the variable to merge by
#' @param addN if TRUE the number of samples will be
#'             added in a column called 'n'
#'
#' @return a phyloseq object with correctly merged
#'         sample data
#' @export
#'
#' @examples
merge_samples_mod <- function(phylo, group, addN = T){
  library(tidyverse)
  
  # Merging using native function
  p <- merge_samples(phylo, group = group, fun = sum)
  
  # Summarising function
  summary_fun <- function(x){
    u <- unique(x)
    if(length(u) > 1){
      vals <- paste0("[", paste(u, collapse = ", "), "]")
      msg <- sprintf("more than one value found: %s, picking the first one.",
                     vals)
      warning(msg, call. = F)
    }
    x[1]
  }
  
  # Summary of all variable, if a variable has more
  # than one value after grouping the firs one will
  # be chosen and a warning will be thrown.
  d <- sample_data(phylo) %>%
    as("data.frame") %>%
    group_by_at(vars(!!group))
  
  if(addN){
    d <- mutate(d, n = n())
  }
  d <- d %>%
    summarise_all(summary_fun) %>%
    data.frame(row.names = sample_names(p))
  
  sample_data(p) <- sample_data(d)
  p
}

getDiversityIndexes <- function(otutable, by = "row", metadata = NULL) {
  # Computes and table the principal biodiversity indexes
  #
  # Args:
  #   otutable: the otutable
  #   by: the position of the OTUs
  #   metadata: a data frame containing
  #             sample information (one for each row)
  #
  # Returns:
  #   A table containing the main biodiversity indexes
  
  library(vegan)
  
  if(by == "row")
    otutable <- data.frame(t(otutable))
  
  N <- rowSums(otutable)
  otus <- specnumber(otutable)
  shannon <- diversity(otutable, index = "shannon")
  invsimpson <- diversity(otutable, index = "invsimpson")
  evenness <- shannon/log(otus)
  sing <- rowSums(round(otutable) == 1)
  doub <- rowSums(round(otutable) == 2)
  chao <- estimateR(round(otutable))["S.chao1",]
  goods <- (1 - (sing/N)) * 100
  
  
  # Tabling everything
  table <- data.frame(Number_of_clones=N,
                      Number_of_OTUs=otus,
                      Number_of_singletons=sing,
                      Number_of_doubletons=doub,
                      Chao1_richness=chao,
                      Shannon_diversity=shannon,
                      InvSimpson=invsimpson,
                      Evenness=evenness,
                      Good_coverage_estimator=goods)
  
  if(!is.null(metadata))
    return(data.frame(table, metadata = metadata))
  return(table)
}

#' PCoA with centroids
#' 
#' Computes a ggplot2 friendly data frame from an ordination.
#'
#' @param ord an ordination. The object must have scores associated
#'             with it, otherwise the function will throw an error.
#' @param meta a data frame. A data frame with metadata.
#' @param sample the name of the column of the data frame specified in meta
#'                with the name of the samples. The name of the samples must
#'                be equal to the names of the rows of the score matrix.
#' @param centroids a character vector with the name of the column of the
#'                   data frame specified in meta to be used for centroids
#'                   calculation. Centroids are calculated using the envfit
#'                   function
#' @param spider if TRUE a dataframe containing spiders will be returned
#' @param hull if TRUE hulls willbe returned
#' @param k a numeric vector used for coordinate selection from ordination
#' @param which a character vector of lenght 1. If the scores function returns
#'               a list (such as with cca or other ecological analyses) this
#'               element will be taken into account.
#' @param ... additional arguments passed to envfit
#'
#' @return a list of data frames
#' @export
#'
#' @examples
#' 
#' library(vegan)
#' library(ggplot2)
#' 
#' x <- matrix(nrow = 100, ncol = 10)
#' colnames(x) <- c(
#'   paste("a", 1:5, sep = "."),
#'   paste("b", 1:5, sep = ".")
#' )
#' for(i in 1:100){
#'   a.values <- rpois(5, lambda = sample(10:1000, 1))
#'   b.values <- rpois(5, lambda = sample(10:1000, 1))
#'   x[i,] <- c(a.values, b.values)
#' }
#' meta <- data.frame(fct = rep(c("a", "b"), each = 5),
#'                    # casual factor
#'                    fct2 = sample(c("D", "E"), ncol(x), replace = T),
#'                    id = colnames(x))
#' x <- prop.table(t(x),1)
#' d <- vegdist(x, "bray")
#' pcoa <- cmdscale(d, eig = T, k = 6)
#' 
#' res <- getOrdination(pcoa, meta, sample = "id", 
#'                      centroids = c("fct", "fct2"), 
#'                      spider = T, hull = T, k = 1:2)
#' 
#' labs <- sprintf("PCoA%d (%.2f%%)", seq_along(res$variance), res$variance)
#' ggplot(res$data, aes(x = x, y = y)) +
#'   # First factor
#'   geom_segment(data = res$spider$fct,
#'                aes(xend = xend, yend = yend, color = fct),
#'                size = .3, alpha = .5, show.legend = F) +
#'   geom_point(aes(color = fct), size = 3) +
#'   geom_label(data = res$centroids$fct,
#'              aes(label = fct, color = fct),
#'              show.legend = F) +
#'   geom_polygon(data = res$hull$fct,
#'                aes(color = fct),
#'                size = .4, linetype = 3,
#'                fill = "transparent") +
#'   # Second factor
#'   geom_segment(data = res$spider$fct2,
#'                aes(xend = xend, yend = yend, color = fct2),
#'                size = .3, alpha = .5, show.legend = F) +
#'   geom_point(aes(color = fct2)) +
#'   geom_label(data = res$centroids$fct2,
#'              aes(label = fct2, color = fct2),
#'              show.legend = F) +
#'   geom_polygon(data = res$hull$fct2,
#'                aes(color = fct2),
#'                size = .2, linetype = 4,
#'                fill = "transparent") +
#'   labs(x = labs[1], y = labs[2])
getOrdination <- function(ord, 
                          meta = NULL, 
                          sample = NULL, 
                          centroids = NULL, 
                          spider = F,
                          hull = F,
                          k = 1:2,
                          which = "sites",
                          ...){
  
  s <- scores(ord)
  if(is.list(s)){
    s <- get(which, s)
  }
  
  data <- data.frame(s[,k])
  colnames(data) <- c("x", "y")
  
  eig <- eigenvals(ord)
  if(!is.null(eig)){
    # Varince explained by axes
    variance <- (eig / sum(eig)) * 100
    variance <- variance[k]
    names(variance) <- colnames(data)
  }
  
  # Formatting data
  data[,sample] <- rownames(data)
  rownames(data) <- NULL
  
  # Merging with metadata
  if(!is.null(meta) & !is.null(sample)){
    data <- merge(data, meta, by = sample)
  }
  
  res <- list(data = data, variance = variance)
  
  # If a group is present, calculate centroids
  if(!is.null(centroids) & !is.null(meta)){
    meta.sub <- subset(meta, select = centroids)
    for(i in ncol(meta.sub)){
      meta.sub[,i] <- as.character(meta.sub[,i])
    }
    
    # Calculate centroids
    fit <- envfit(ord, meta.sub, choices = k, ...)
    centroids.data <- data.frame(fit$factors$centroids)
    colnames(centroids.data) <- names(variance)
    
    # Formatting datasets
    var.names <- colnames(meta.sub)
    centroids.data <- lapply(var.names, function(v){
      r.names <- paste0(v, unique(get(v, data)))
      sub.centr <- centroids.data[r.names,]
      rgx <- paste0("^", v)
      sub.centr[,v] <- gsub(rgx, "", rownames(sub.centr))
      rownames(sub.centr) <- NULL
      sub.centr
    })
    names(centroids.data) <- var.names
    res$centroids <- centroids.data
    
    if(spider){
      # Adding centroids coordinates for 
      # spider plotting
      spiders <- lapply(centroids.data, function(c){
        v.names <- colnames(c)
        sp <- subset(data, select = v.names)
        colnames(c)[1:2] <- c("xend", "yend")
        merge(sp, c, by = v.names[3])
      })
      res$spider <- spiders
    }
    
    if(hull){
      # Adding hull coordinates for hulls
      hulls <- lapply(centroids, function(c){
        h <- subset(data, select = c(c, names(variance)))
        grp <- get(c, h)
        grp.h <- lapply(split(h, grp), function(x){
          i <- chull(subset(x, select = names(variance)))
          x[i,]
        })
        do.call(rbind, grp.h)
      })
      names(hulls) <- centroids
      res$hull <- hulls
    }
  }
  return(res)
}

pairwise.adonis <- function(x,factors, sim.method, p.adjust.m, ...) {
  # Computes and table post-hoc test for adonis
  #
  # Args:
  #   x: the otutable
  #   factors: factors to test differences
  #   sim.method: similarity index (to be passed to adonis)
  #   p.adjust.m: pvalue adjust method (to be passed to p.adjust)
  #
  # Returns:
  #   A table containing post-hoc tests
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] , 
                method =sim.method, ...)
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}
