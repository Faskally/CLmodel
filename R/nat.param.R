# copy of non exported function from mgcv
nat.param <- function (X, S, rank = NULL, type = 0, tol = .Machine$double.eps^0.8,
    unit.fnorm = TRUE)
{
    if (type == 2 || type == 3) {
        er <- eigen(S, symmetric = TRUE)
        if (is.null(rank) || rank < 1 || rank > ncol(S)) {
            rank <- sum(er$value > max(er$value) * tol)
        }
        null.exists <- rank < ncol(X)
        E <- rep(1, ncol(X))
        E[1:rank] <- sqrt(er$value[1:rank])
        X <- X %*% er$vectors
        col.norm <- colSums(X^2)
        col.norm <- col.norm/E^2
        av.norm <- mean(col.norm[1:rank])
        if (null.exists)
            for (i in (rank + 1):ncol(X)) {
                E[i] <- sqrt(col.norm[i]/av.norm)
            }
        P <- t(t(er$vectors)/E)
        X <- t(t(X)/E)
        if (null.exists && type == 3 && rank < ncol(X) - 1) {
            ind <- (rank + 1):ncol(X)
            rind <- ncol(X):(rank + 1)
            Xn <- X[, ind, drop = FALSE]
            n <- nrow(Xn)
            one <- rep(1, n)
            Xn <- Xn - one %*% (t(one) %*% Xn)/n
            um <- eigen(t(Xn) %*% Xn, symmetric = TRUE)
            X[, rind] <- X[, ind, drop = FALSE] %*% um$vectors
            P[, rind] <- P[, ind, drop = FALSE] %*% (um$vectors)
        }
        if (unit.fnorm) {
            ind <- 1:rank
            scale <- 1/sqrt(mean(X[, ind]^2))
            X[, ind] <- X[, ind] * scale
            P[ind, ] <- P[ind, ] * scale
            if (null.exists) {
                ind <- (rank + 1):ncol(X)
                scalef <- 1/sqrt(mean(X[, ind]^2))
                X[, ind] <- X[, ind] * scalef
                P[ind, ] <- P[ind, ] * scalef
            }
        }
        else scale <- 1
        return(list(X = X, D = rep(scale^2, rank), P = P, rank = rank,
            type = type))
    }
    qrx <- qr(X, tol = .Machine$double.eps^0.8)
    R <- qr.R(qrx)
    RSR <- forwardsolve(t(R), t(forwardsolve(t(R), t(S))))
    er <- eigen(RSR, symmetric = TRUE)
    if (is.null(rank) || rank < 1 || rank > ncol(S)) {
        rank <- sum(er$value > max(er$value) * tol)
    }
    null.exists <- rank < ncol(X)
    D <- er$values[1:rank]
    X <- qr.Q(qrx, complete = FALSE) %*% er$vectors
    P <- backsolve(R, er$vectors)
    if (type == 1) {
        E <- c(sqrt(D), rep(1, ncol(X) - length(D)))
        P <- t(t(P)/E)
        X <- t(t(X)/E)
        D <- D * 0 + 1
    }
    if (unit.fnorm) {
        ind <- 1:rank
        scale <- 1/sqrt(mean(X[, ind]^2))
        X[, ind] <- X[, ind] * scale
        P[, ind] <- P[, ind] * scale
        D <- D * scale^2
        if (null.exists) {
            ind <- (rank + 1):ncol(X)
            scalef <- 1/sqrt(mean(X[, ind]^2))
            X[, ind] <- X[, ind] * scalef
            P[, ind] <- P[, ind] * scalef
        }
    }
    list(X = X, D = D, P = P, rank = rank, type = type)
}
