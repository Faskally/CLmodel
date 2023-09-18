
#' Estimate capture probabilites from electrofishing data
#'
#' This function uses the marginal likelihood of capture probabilities
#' to estimate model parameters
#'
#' @param object a gmrf smooth object
#' @param data a data.frame containing all relavent info
#' @param knots todo
#'
#' @return glm type object
#'
#' @export
smooth.construct.gmrf.smooth.spec <- function(object, data, knots) {
    k <- factor(rownames(object$xt$penalty), levels = rownames(object$xt$penalty))
    x <- data[[object$term]]
    x <- factor(x, levels = levels(k))
    #k <- factor(levels(x), levels = levels(x))
    if (object$bs.dim < 0) object$bs.dim <- length(levels(k))
    if (object$bs.dim > length(levels(k))) stop("GMRF basis dimension set too high")
    object$X <- model.matrix(~x - 1, )
    # penalty
    object$S[[1]] <- object$xt$penalty
    if (ncol(object$S[[1]]) != nrow(object$S[[1]]))
        stop("supplied penalty not square!")
    if (ncol(object$S[[1]]) != ncol(object$X))
        stop("supplied penalty wrong dimension!")
    if (!is.null(colnames(object$S[[1]]))) {
        a.name <- colnames(object$S[[1]])
        if (all.equal(levels(k), a.name) != TRUE) {
            stop("penalty column names don't match supplied area names!")
        }
        else {
            if (all.equal(sort(a.name), a.name) != TRUE) {
              object$S[[1]] <- object$S[[1]][levels(k), ]
              object$S[[1]] <- object$S[[1]][, levels(k)]
            }
        }
    }

    if (object$bs.dim < length(levels(k))) {
        mi <- which(colSums(object$X) == 0)
        np <- ncol(object$X)
        if (length(mi) > 0) {
            object$X <- rbind(matrix(0, length(mi), np), object$X)
            for (i in 1:length(mi)) object$X[i, mi[i]] <- 1
        }
        rp <- nat.param(object$X, object$S[[1]], type = 0)
        ind <- (np - object$bs.dim + 1):np
        object$X <- if (length(mi))
            rp$X[-(1:length(mi)), ind]
        else rp$X[, ind]
        object$P <- rp$P[, ind]
        object$S[[1]] <- diag(c(rp$D[ind[ind <= rp$rank]], rep(0,
            sum(ind > rp$rank))))
        object$rank <- rp$rank
    }
    else if (!is.null(object $ xt $ rank)) {
      object$rank <- object $ xt $ rank
    }
    else {
        ev <- eigen(object$S[[1]], symmetric = TRUE, only.values = TRUE)$values
        object$rank <- sum(ev > .Machine$double.eps^0.8 * max(ev))
    }
    object$null.space.dim <- ncol(object$X) - object$rank
    object$knots <- k
    object$df <- ncol(object$X)
    object$plot.me <- FALSE
    #object$fixed <- FALSE # force penalty - no sense in allowing fixed to be false
    # specify the constraints
    #object $ C <- NULL

    class(object) <- "gmrf.smooth"
    object
}


#' Estimate capture probabilites from electrofishing data
#'
#' This function uses the marginal likelihood of capture probabilities
#' to estimate model parameters
#'
#'
#' @param object a gmrf smooth object
#' @param data a data.frame containing all relavent info
#'
#' @return glm type object
#'
#' @export
Predict.matrix.gmrf.smooth <- function (object, data) {
    x <- factor(data[[object$term]], levels = levels(object$knots))
    X <- model.matrix(~x - 1)
    if (!is.null(object$P))
        X <- X %*% object$P
    X
}
