#' @import methods
NULL

#' DataSet for hypergeometric post-hoc test
#'
#' Main class definition.
#'
#' @slot data numeric. A matrix with abundance data
#' @slot fct character. A character vector with categories
#'
#' @return A DataSet object
#' @export
#'
#' @examples
#' x <- matrix(rnorm(1000), ncol = 10)
#' cat <- sample(letters[1:4], 10, replace = T)
#'
#' data <- new("DataSet", data = x, fct = cat)
DataSet <- setClass("DataSet",
                    slots = c(
                      data = "matrix",
                      fct = "character"
                    ))

#' Validation funciton for DataSet object
#'
#' @param object the [DataSet] object
#'
#' @return TRUE if parameters are correct
#' @export
validDataSetObject <- function(object){
  if(!is.numeric(object@data)){
    "data must be numeric"
  }else if(ncol(object@data) != length(object@fct)){
    "fct != ncol(data)"
  }else{
    TRUE
  }
}
setValidity("DataSet", validDataSetObject)

# Setting old S3 class
setOldClass("hclust")

# Setting hclust or NULL
setClassUnion("hclustOrNULL",members=c("hclust", "NULL"))

#' A grouped DataSet
#'
#' This class must be constructed from a \code{\link{DataSet}} object
#' using the \code{\link{groupData}} funciton
#'
#' @slot tree hclust.
#' @slot grp numeric.
#' @slot grp.method character. The grouping method used
#'
#' @return a GroupedDataSet object
#'
#' @examples
#' x <- data.frame(matrix(rnorm(1000), ncol = 10))
#' cat <- sample(letters[1:4], 10, replace = T)
#'
#' data <- buildDataSet(x, cat)
#'
#' grouped <- groupData(data,
#'                      dist.method = "manhattan",
#'                      clust.method = "complete",
#'                      cut.method = "cutreeHybrid")
#' class(grouped)
GroupedDataSet <- setClass("GroupedDataSet", contains = "DataSet",
                           slots = c(
                             tree = "hclustOrNULL",
                             grp = "numeric",
                             grp.method = "character"
                           ))

#' Validation funciton for GroupedDataSet object
#'
#' @param object the [GroupedDataSet] object
#'
#' @return TRUE if parameters are correct
#' @export
validGroupedDataSetObject <- function(object){
  if(length(object@grp) != nrow(object@data)){
    "grp != nrow(data)"
  }else{
    TRUE
  }
}
setValidity("GroupedDataSet", validGroupedDataSetObject)

#' Constructor for DataSet objects
#'
#' Main constructor for [DataSet] objects. If the \code{grp} parameter
#' is specified an object a [GroupedDataSet] is returned. Each group
#' must be reported with a number whereas elements not assigned to
#' any group must be reported with \code{NA}.
#'
#' @param data numeric. A matrix with abundance data
#' @param fct character. A character vector with categories
#' @param grp character. A numeric vector with groups. If \code{NULL}
#' no group will be assigned
#'
#' @return A DataSet object
#' @export
#'
#' @examples
#' x <- matrix(rnorm(1000), ncol = 10)
#' cat <- sample(letters[1:4], 10, replace = T)
#'
#' data <- buildDataSet(data = x, fct = cat)
#'
#' x <- as.data.frame(x)
#' cat <- as.factor(cat)
#' data <- buildDataSet(data = x, fct = cat)
buildDataSet <- function(data, fct, grp = NULL){
  if(!is.matrix(data)){
    message("converting data to matrix")
    data <- as.matrix(data)
  }
  
  if(!is.character(fct)){
    message("converting fct to character")
    fct <- as.character(fct)
  }
  
  if(!is.null(grp) & !is.numeric(grp)){
    message("converting grp to character")
    grp <- as.numeric(grp)
  }
  
  res <- DataSet(data = data, fct = fct)
  if(is.null(grp)){
    res
  }else{
    GroupedDataSet(res, grp = grp, tree = NULL)
  }
}


#' Title
#'
#' @param obj
#' @param categories
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
setGeneric("hpostSingle", function(obj, categories) standardGeneric("hpostSingle"))

setMethod("hpostSingle", "GroupedDataSet", function(obj, categories){
  stopifnot(length(categories) == nrow(obj@data))
  
  pop.size <- nrow(obj@data)
  
  cat.all <- table(categories) %>%
    dplyr::as_tibble()
  
  dplyr::tibble(grp = obj@grp,
                categories = categories) %>%
    na.omit() %>%
    dplyr::group_by(grp, categories) %>%
    dplyr::summarise(q = dplyr::n()) %>%
    dplyr::mutate(k = sum(q)) %>%
    dplyr::left_join(cat.all, by = "categories") %>%
    dplyr::mutate(m = n,
                  n = pop.size - m,
                  fraction.pop = m / (n + m),
                  fraction.obs =  q / k,
                  exp.q = k * fraction.pop,
                  log2.fold.enrichment = log2(fraction.obs / fraction.pop),
                  p.value = ifelse(log2.fold.enrichment > 0,
                                   phyper(q = q - 1, m = m, n = n, k = k, lower.tail = F),
                                   phyper(q = q, m = m, n = n, k = k, lower.tail = T)),
                  adj.p.value = p.adjust(p.value, "BH"))
})
