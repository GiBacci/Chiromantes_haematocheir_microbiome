#' Prepare data for sunburs plotting
#' 
#' This function takes a hierarchical matrix as input
#' and transform it into a data frame that can be used
#' in combination with geom_rect and coord_polar of 
#' ggplot2 package.
#' 
#' The columns fo the data frame are treated as hierarchical
#' levels, from the first column (the highest level) to the
#' last (the most specific level)
#' 
#'
#' @param data a data.frame
#' @param bar.width the width of the bars
#' @param collapse if TRUE consecutive layers with the same annotation
#'                 will be collapsed togheter (if used with a bar.width < 1
#'                 results can be graphically clumsy)
#'
#' @return a data frame that can be used in ggplot2
#' @export
#'
#' @examples
#' dat <- data.frame(level1 = rep(c("a", "b"), each=3),
#'                   level2 = paste0(rep(c("a", "b"), each=3), 1:3),
#'                   stringsAsFactors = FALSE)
#' 
#' res <- getSunburstData(dat)
#' ggplot(res, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
#'   geom_rect(color = "white", show.legend = F) +
#'   coord_polar("y") +
#'   xlim(-1, NA)
getSunburstData <- function(data, bar.width = 1, collapse = F){
  
  if(bar.width > 1 | bar.width < 0){
    stop("bar.width must be comprised between 0 and 1")
  }
  
  # Get last elements (useful for labelling)
  lasts <- unique(data)
  lasts <- apply(lasts, 1, function(x){
    no.na <- na.omit(x)
    last <- no.na[length(no.na)]
    
    found <- sum(lasts[[length(no.na)]] == last, na.rm = T)
    if(found > 1)
      return(NA)
    last
  })
  lasts <- na.omit(lasts)
  
  # Transforming NA into explicit levels
  for(i in 2:ncol(data)){
    x <- data[,i]
    y <- data[,i-1]
    extr <- is.na(x)
    x[extr] <- paste0("NA__", y[extr])
    data[,i] <- x
  }
  
  # Ordering based on occurrence
  for(i in 1:ncol(data)){
    lvl <- names(sort(table(data[[i]]), decreasing = T))
    data[[i]] <- factor(data[[i]], levels = lvl)
  }
  
  # Re-ordering based on all columns (from first to last)
  data <- data[do.call(order, data),]
  
  # Building dataset for plotting
  l <- lapply(1:ncol(data), function(i){
    x <- data[[i]]
    x <- table(x)[unique(x)]
    
    to <- cumsum(x)
    from <- c(0, to[-length(to)])
    data.frame(label = names(to), 
               ymin = from, 
               ymax = to, 
               xmin = i,
               xmax = i + bar.width,
               stringsAsFactors = F, 
               row.names = NULL)
  })
  res <- do.call(rbind, l)
  
  if(collapse){
    if(bar.width < 1){
      warning("Collapsing same tax levels with bar.width < 1")
    }
    res <- split(res, res$label)
    res <- lapply(res, function(x){
      if(nrow(x) > 1){
        x$xmin <- min(x$xmin)
        x$xmax <- max(x$xmax)
        return(x[1,])
      }
      x
    })
    res <- do.call(rbind, res)
  }
  res <- data.frame(res, row.names = NULL,
                    stringsAsFactors = F)
  
  # Removing NAs
  nas <- grep("NA__", res$label)
  if(length(nas) > 0) res <- res[-nas,]
  
  # Determing x range for label angles
  x.range <- c(min(res$ymin), max(res$ymax))
  x.length <- x.range[2] - x.range[1]
  
  # Adding y and x values for labels ( (ymin - ymax) / 2 )
  res <- transform(res,
                   lab.y = ymin + ((ymax - ymin) / 2),
                   lab.x = xmin + ((xmax - xmin) / 2),
                   last = label %in% lasts)
  
  # Adding angle for labels
  transform(res,
            angle = (lab.y/x.length)*360)
}


#' Get sunburst labels
#' 
#' This function return a geom_text layer with
#' labels formatted according to poola coordinates.
#' 
#' Data must contain:
#'  1. A columns with labels
#'  2. A columns with angles (normally computed as a 
#'     fraction of 360 based on x axis, e.g. (x/max(x)) * 360)
#'  3. A column with x position of labels
#'  4. A column with y position of labels
#'
#' @param data a data frame with columns described above
#' @param alignment the alignment of the labels in respect 
#'                  with poolar coordinates
#' @param position the position of the labels
#' @param offset the offset (only if position == outside)
#' @param angle the name of the column with angles
#' @param xmax the name of the column with x values
#' @param label the name of the column with labels
#' @param lab.y the name of the column with y value
#' @param condensed if TRUE condensed label will be used
#' @param condensed.legend if TRUE a legend with condensed label is added
#' @param ... additional options passed ot geom_text
#'
#' @return a ggplot2 layer with labels
#' @export
#'
#' @examples
#' dat <- data.frame(level1 = rep(c("a", "b"), each=3),
#'                   level2 = paste0(rep(c("a", "b"), each=3), 1:3),
#'                   stringsAsFactors = FALSE) 
#'                   
#' res <- getSunburstData(dat)
#' 
#' ggplot(res, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
#'   geom_rect(color = "white", show.legend = F) +
#'   coord_polar("y") +
#'   xlim(-1, NA) +
#'   polar_labels(subset(res, last == T), alignment = "perpendicular",
#'                position = "outside") +
#'   polar_labels(subset(res, last == F), alignment = "normal",
#'                position = "inside")
polar_labels <- function(data, 
                         alignment = c("perpendicular", "tangent", "normal"), 
                         position = c("outside", "inside"), 
                         offset = .1,
                         angle = "angle",
                         xmax = "xmax",
                         label = "label",
                         lab.y = "lab.y",
                         condensed = F,
                         condensed.legend = F,
                         use.shadowtext = F,
                         ...){
  
  # get alignment type and position
  alignment <- alignment[1]
  position <- position[1]
  
  # get data
  angle.ind <- which(names(data) == angle)
  xmax.ind <- which(names(data) == xmax)
  
  angle <- data[[angle.ind]]
  xmax <- data[[xmax.ind]]
  
  # if alignment is perpendicular labels
  # are rotated and flipped. If is tangent labels
  # are partially rotated and flipped.
  if(alignment == "perpendicular"){
    angle <- 90 - angle
    hjust <- ifelse( angle < -90, 1, 0)
    angle <- ifelse( angle < -90, angle + 180, angle)
  }else if(alignment == "tangent"){
    angle <- - angle
    angle <- ifelse( angle < -90 & angle > -270, angle + 180, angle)
    hjust <- .5
  }else if(alignment == "normal"){
    angle <- rep(0, length(angle))
    hjust <- NA
  }else{
    stop("please specify a valid alignment option")
  }
  
  # Adding angle to data
  data[[angle.ind]] <- angle
  data <- transform(data, hjust = hjust)
  
  # Add position
  if(position == "outside"){
    data <- transform(data, lab.x = xmax + offset)
  }else if(position == "inside"){
    data <- data
  }else{
    stop("please specify a valid position option")
  }
  
  
  if(condensed){
    index <- ceiling(1:nrow(data)/length(letters))
    
    l <- if(nrow(data) > length(letters)){
      letters
    }else{
      letters[1:nrow(data)]
    }
    
    data <- transform(data, condensed = paste0(l, index))
    
    labels <- get(label, data)
    label <- "condensed"
    breaks <- get(label, data)
  }
  
  res <- if(use.shadowtext){
    shadowtext::geom_shadowtext(data = data, 
              aes(x = lab.x, 
                  label = !!ensym(label), 
                  y = !!ensym(lab.y), 
                  angle = angle, 
                  hjust = hjust),
              ...)
  }else{
    geom_text(data = data, 
              aes(x = lab.x, 
                  label = !!ensym(label), 
                  y = !!ensym(lab.y), 
                  angle = angle, 
                  hjust = hjust),
              ...)
  }
  
  if(condensed.legend){
    return(
      list(
        geom = res,
        legend = scale_discrete_identity(
          aesthetics = "label",
          breaks = breaks,
          labels = labels,
          guide = "legend"
        )
      )
    )
  }else{
    return(res)
  }
}


#' Closure for returning a key_glyph formatting function
#' 
#' This function will return a key_glyph formatting function based
#' on the text sise provided
#'
#' @param lab.size the size of the text to be included in the legend
#'
#' @return a key_glyph function to be included in the geom_text
#' @export
#'
#' @examples
#' dat <- data.frame(level1 = rep(c("a", "b"), each=3),
#'                   level2 = paste0(rep(c("a", "b"), each=3), 1:3),
#'                   stringsAsFactors = FALSE) 
#'                   
#' res <- getSunburstData(dat)
#' 
#' key_glyph <- get_draw_labels(8)
#' ggplot(res, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
#'   geom_rect(color = "white", show.legend = F) +
#'   coord_polar("y") +
#'   xlim(-1, NA) +
#'   polar_labels(res, alignment = "normal", position = "inside", 
#'                condensed = T, condensed.legend = T,
#'                key_glyph = key_glyph) +
#'   theme_void(base_size = 8)
get_draw_labels <- function(lab.size){
  
  function(data, params, size) {
    grid::textGrob(
      paste0(data$label, ":"), 1., 0.5,
      hjust = 1, # right justified
      rot = data$angle,
      gp = grid::gpar(
        col = scales::alpha(data$colour, data$alpha),
        fontfamily = data$family,
        fontface = data$fontface,
        fontsize = .8*lab.size # match font size to label font
      )
    )
  }
}