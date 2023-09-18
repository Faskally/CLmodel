#' Modelling Framework for the Estimation of Salmonid Abundance in Scottish Rivers
#'
#' @docType package
#' @name CLmodel
#' @description Functions to estimate juvenile densities of Salmon and Trout on river networks.
#' @import mgcv
#'
#' @importFrom grDevices grey
#' @importFrom graphics lines locator points title
#' @importFrom stats AIC BIC binomial model.matrix rnorm
#' @importFrom utils capture.output tail globalVariables
#' @importFrom igraph E V graph.edgelist is.connected degree delete.vertices get.shortest.paths V<- E<-
#' @importFrom ggplot2 theme element_blank geom_line aes autoplot
#' @importFrom ggmap ggmap get_map
#' @importFrom rgeos gIntersection
#' @importFrom sp spTransform SpatialPoints CRS
#' @importFrom FNN get.knnx
#' @importFrom OpenStreetMap openmap
#'
NULL