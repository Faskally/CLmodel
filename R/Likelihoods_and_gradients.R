
#
# Purpose: fit capture probabilities and abundances to electrofishing data
# 
# author: CP Millar, millarc@marlab.ac.uk
# origional date: 12/2014
#

#' Estimate capture probabilites from electrofishing data
#'
#' This function uses the marginal likelihood of capture probabilities
#' to estimate model parameters 
#' 
#'
#' @param formula a formula object
#' @param data a data.frame containing all relavent info
#' @return glm type object
#' @export
#' @examples
#' # none yet
# efp_old <- function(formula, lformula = NULL, data = NULL, usegrad = TRUE, hessian = FALSE) {

#   if (is.null(data)) stop("must suply data")
#   data0 <- subset(data, T > 0)

#   # set up model
#   if (nrow(data0) == 1) {
#     G <- matrix(1, 1, 1)
#   } else {
#     Gsetup <- gam(formula, data = data0, fit = FALSE)
#     G <- Gsetup $ X
#   }

#   # objective
#   fun <- function(par, G, T, X, s) {
#     Gb <- G %*% par
#     eGb <- exp(Gb)
#     R <- (s - 1 - X/T)
  
#     loglik_i <- T * (Gb - (R+1) * log(1 + eGb) - log(1 - (1 + eGb)^(-s)) )
#     -1 * sum(loglik_i)
#   }

#   # gradient
#   dfun <- function(par, G, T, X, s) {
#     eGb <- exp(t(t(G) * par))
#     peGb <- 1 + eGb 

#     -1 * colSums( 
#         T * (G / peGb  - 
#           (s - 1 - X/T) * G * eGb / peGb -
#             s * G * eGb / (peGb^(s+1) - peGb)
#         ))
#   }

  
#   par0 <- rep(0, ncol(G))
#   #suppressWarnings(
#   #  opt <- optim(rep(0, ncol(G)), fun, G = G, 
#   #             T = data $ T, X = data $ X, s = data $ s,
#   #             method = "BFGS"))

#   if (usegrad) {
#     opt <- nlminb(start = par0, objective = fun, gradient = dfun,
#                   T = data0 $ T, X = data0 $ X, s = data0 $ s, G = unname(G))
#     opt $ grad <- dfun(opt $ par, unname(G), data0 $ T, data0 $ X, data0 $ s)
#   } else {
#     opt <- nlminb(start = par0, objective = fun,
#                   T = data0 $ T, X = data0 $ X, s = data0 $ s, G = unname(G))    
#   }

#   # predict p for rest of data: NB problems can arise here
#   data $ p <- transpar(opt $ par, gam(formula, data = data, fit = FALSE) $ X)

#   # get loglik conditional on a lambda model
#   if (is.null(lformula)) { # then fit full model

#     # this is the concentrated likelihood
#     loglik <- function(p, T, X, S) {
#       sum(ifelse(T > 0,
#         T*log(T) + T*log(p) + (T*(S - 1) - X)*log(1-p) - T*log(1-(1-p)^S) - T,
#         0))
#     }
#     lmodel <- T ~ 1
#     llik <- loglik(data $ p, data $ T, data $ X, data $ s)
 
#     hessFunc <- function(par) {
#       p <- transpar(par, gam(formula, data = data, fit = FALSE) $ X)
#       loglik(p, data $ T, data $ X, data $ s)  
#     }

#     if (hessian) {
#       H <- numDeriv::hessian(hessFunc, opt $ par)
#       Vb <- solve(-1 * H)
#     }

#   } else {

#     stop("models for lambda not implemented yet in efp()")

#     # get lambdas: they have a gamma distribution (overly complicated I know!)
#     apar <- data $ T + 1
#     bpar <- with(data, (1 - (1-p)^Runs)) # dont care about area for now, only likelihood
#     #data $ density <- apar / bpar # mean ... never on boundary ...
#     data $ density <- (apar-1) / bpar # mode - i.e. ML estimate, can be zero
#     #data $ density <- with(data, T / ((1 - (1-p)^Runs))) # ML estimator - i.e. ML estimate, can be zero
    
#     # this is the full likelihood
#     loglik <- function(p, lambda, T, X, S) {
#       sum(ifelse(lambda > 0, # because if lambda == 0, all T informing the esitmate are 0
#         T*log(p) + (T*(S - 1) - X)*log(1-p) + T*log(lambda) - (1-(1-p)^S) * lambda,
#         0))
#     }
#     llik <- loglik(data $ p, data $ density, data $ T, data $ X, data $ s)
#   }
  
#   opt $ formula <- formula # for printing and summary
#   opt $ lformula <- lformula
#   opt $ llik <- llik
#   opt $ y <- data0 $ X
#   opt $ terms <- Gsetup $ terms
#   opt $ call <- match.call()
#   opt $ aic <- -2 * opt $ llik + 2 * ncol(G)
#   opt $ G <- G
#   opt $ call <- match.call()
#   opt $ coefficients <- opt $ par
#   names(opt $ coefficients) <- colnames(G)
#   opt $ df.null <- nrow(G)
#   opt $ df.residual <- nrow(G) - ncol(G)
#   opt $ rank <- ncol(G)
#   opt $ fitted <- data $ p
#   opt $ residuals <- rep(0, nrow(data0))
#   opt $ null.deviance <- NA
#   opt $ deviance <- NA 
#   opt $ family <- binomial()
#   if (hessian) {
#     opt $ Vb <- Vb
#   }
#   class(opt) <- c("efp", "glm", "lm")
#   opt
# }



#' Estimate capture probabilites from electrofishing data
#'
#' This function uses the marginal likelihood of capture probabilities
#' to estimate model parameters 
#' 
#'
#' @param formula a formula object
#' @param data a data.frame containing all relavent info
#' @return glm type object
#' @export
#' @examples
#' # none yet
efp <- function(formula, data = NULL, passes = NULL, verbose=TRUE, init = "0") {

  if (!exists("stanmod")) {
    message("Building optimiser for first use...")
    stanmod <- rstan::stan_model(model_code = "
      data {
        int<lower=0> N; // number of observations
        int<lower=0> K; // number of parameters
        real S[N]; // the number of fishing passes
        real R[N]; // Zippins R (see seber p 312, eq 7.22)
        real T[N]; // total catches
        matrix[N,K] A;
      }
      parameters {
        vector[K] alpha;
      } 
      model {
        vector[N] expeta;
        expeta <- exp(A * alpha);
        for (i in 1:N) {
          real p;
          p <- expeta[i]/(1.0 + expeta[i]);
          increment_log_prob(T[i] * log(p));
          increment_log_prob(T[i] * R[i] * log(1-p));
          increment_log_prob(-T[i] * log(1 - (1-p)^S[i]) );
        }
      }")
    assign("stanmod", stanmod, .GlobalEnv)
  }

  if (is.null(data)) stop("must supply data")
  data0 <- subset(data, T > 0)

  if (is.null(passes)) stop("must supply the number of fishing runs")
  data0 $ S <- data0[[passes]]

  # set up model
  if (nrow(data0) == 1) {
    G <- matrix(1, 1, 1)
  } else {
    Gsetup <- gam(formula, data = data0, fit = FALSE)
    G <- Gsetup $ X
  }

  standat <- 
    list(N = nrow(G), K = ncol(G), 
         S = data0 $ S, T = data0 $ T, R = with(data0, S - 1 - Z),
         A = G)

  opt <- optimizing(stanmod, data = standat, algorith = "BFGS", hessian = TRUE, verbose = verbose, init = init)

  # predict p for rest of data: NB problems can arise here
  data $ p <- try(transpar(opt $ par, gam(formula, data = data, fit = FALSE) $ X))
   
  opt $ formula <- formula # for printing and summary
  opt $ llik <- opt $ value
  opt $ terms <- Gsetup $ terms
  opt $ call <- match.call()
  opt $ aic <- -2 * opt $ llik + 2 * ncol(G)
  opt $ G <- G
  opt $ coefficients <- opt $ par
  names(opt $ coefficients) <- colnames(G)
  opt $ df.null <- nrow(G)
  opt $ df.residual <- nrow(G) - ncol(G)
  opt $ rank <- ncol(G)
  opt $ fitted <- data $ p
  opt $ residuals <- rep(0, nrow(data0))
  opt $ null.deviance <- NA
  opt $ deviance <- NA 
  opt $ family <- binomial()
  opt $ Vb <- solve(-1 * opt $ hessian)
  opt $ Gsetup <- Gsetup

  # get a gam container
  # g1 <- gam(G = Gsetup)
  # g1 $ coefficients[] <- opt $ par       
  # g1 $ Vp[] <- opt $ Vb
  # g1 $ family <- binomial()
  # X <- predict(g1, type = "lpmatrix")
  # g1 $ linear.predictors <-  c(X %*% g1 $ coef)
  # g1 $ fitted.values <- c(1/(1 + exp(-g1 $ linear.predictors)))
  # g1 $ aic <- opt $ aic

  class(opt) <- c("efp", "glm", "lm")
  opt
}


#' Estimate capture probabilites from electrofishing data
#'
#' This function uses the marginal likelihood of capture probabilities
#' to estimate model parameters 
#' 
#'
#' @param formula a formula object
#' @param data a data.frame containing all relavent info
#' @return glm type object
#' @export
#' @examples
#' # none yet
efp.fit <- function(X, y, data, ...) {

}




#' Estimate capture probabilites from electrofishing data
#'
#' This function uses the marginal likelihood of capture probabilities
#' to estimate model parameters 
#' 
#'
#' @param formula a formula object
#' @param data a data.frame containing all relavent info
#' @return glm type object
#' @export
#' @examples
#' # none yet
addterm.efp <- function (object, scope, scale = 0, test = c("none", "Chisq", 
    "F"), k = 2, sorted = FALSE, trace = FALSE, ...) 
{
    Fstat <- function(table, rdf) {
        dev <- table$Deviance
        df <- table$Df
        diff <- pmax(0, (dev[1L] - dev)/df)
        Fs <- diff/(dev/(rdf - df))
        Fs[df < .Machine$double.eps] <- NA
        P <- Fs
        nnas <- !is.na(Fs)
        P[nnas] <- safe_pf(Fs[nnas], df[nnas], rdf - df[nnas], 
            lower.tail = FALSE)
        list(Fs = Fs, P = P)
    }
    if (missing(scope) || is.null(scope)) 
        stop("no terms in scope")
    if (!is.character(scope)) 
        scope <- add.scope(object, update.formula(object, scope))
    if (!length(scope)) 
        stop("no terms in scope for adding to object")
    oTerms <- attr(terms(object), "term.labels")
    int <- attr(object$terms, "intercept")
    ns <- length(scope)
    dfs <- dev <- numeric(ns + 1)
    names(dfs) <- names(dev) <- c("<none>", scope)
    add.rhs <- paste(scope, collapse = "+")
    add.rhs <- eval(parse(text = paste("~ . +", add.rhs)))
    new.form <- update.formula(object, add.rhs)
    oc <- object$call
    Terms <- terms(new.form)
    oc$formula <- Terms
    fob <- list(call = oc, terms = Terms)
    class(fob) <- class(object)
    x <- model.matrix(Terms, model.frame(fob, xlev = object$xlevels), 
        contrasts = object$contrasts)
    n <- nrow(x)
    oldn <- length(object$residuals)
    y <- object$y
    newn <- length(y)
    if (newn < oldn) 
        warning(sprintf(ngettext(newn, "using the %d/%d row from a combined fit", 
            "using the %d/%d rows from a combined fit"), newn, 
            oldn), domain = NA)
    wt <- object$prior.weights
    if (is.null(wt)) 
        wt <- rep(1, n)
    Terms <- attr(Terms, "term.labels")
    asgn <- attr(x, "assign")
    ousex <- match(asgn, match(oTerms, Terms), 0L) > 0L
    if (int) 
        ousex[1L] <- TRUE
    X <- x[, ousex, drop = FALSE]
    z <- efp.fit(X, y, object $ data)
    dfs[1L] <- z$rank
    dev[1L] <- z$deviance
    sTerms <- sapply(strsplit(Terms, ":", fixed = TRUE), function(x) paste(sort(x), 
        collapse = ":"))
    for (tt in scope) {
        if (trace) {
            message(gettextf("trying + %s", tt), domain = NA)
            utils::flush.console()
        }
        stt <- paste(sort(strsplit(tt, ":")[[1L]]), collapse = ":")
        usex <- match(asgn, match(stt, sTerms), 0L) > 0L
        X <- x[, usex | ousex, drop = FALSE]
        z <- efp.fit(X, y, object $ data)
        dfs[tt] <- z$rank
        dev[tt] <- z$deviance
    }
    if (is.null(scale) || scale == 0) 
        dispersion <- summary(object, dispersion = NULL)$dispersion
    else dispersion <- scale
    fam <- object$family$family
    if (fam == "gaussian") {
        if (scale > 0) 
            loglik <- dev/scale - n
        else loglik <- n * log(dev/n)
    }
    else loglik <- dev/dispersion
    aic <- loglik + k * dfs
    aic <- aic + (extractAIC(object, k = k)[2L] - aic[1L])
    dfs <- dfs - dfs[1L]
    dfs[1L] <- NA
    aod <- data.frame(Df = dfs, Deviance = dev, AIC = aic, row.names = names(dfs), 
        check.names = FALSE)
    o <- if (sorted) 
        order(aod$AIC)
    else seq_along(aod$AIC)
    if (all(is.na(aic))) 
        aod <- aod[, -3]
    test <- match.arg(test)
    if (test == "Chisq") {
        dev <- pmax(0, loglik[1L] - loglik)
        dev[1L] <- NA
        LRT <- if (dispersion == 1) 
            "LRT"
        else "scaled dev."
        aod[, LRT] <- dev
        nas <- !is.na(dev)
        dev[nas] <- safe_pchisq(dev[nas], aod$Df[nas], lower.tail = FALSE)
        aod[, "Pr(Chi)"] <- dev
    }
    else if (test == "F") {
        if (fam == "binomial" || fam == "poisson") 
            warning(gettextf("F test assumes 'quasi%s' family", 
                fam), domain = NA)
        rdf <- object$df.residual
        aod[, c("F value", "Pr(F)")] <- Fstat(aod, rdf)
    }
    aod <- aod[o, ]
    head <- c("Single term additions", "\nModel:", deparse(formula(object)))
    if (scale > 0) 
        head <- c(head, paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}











#' Calculate the relavent statistics for fitting electrifishing models
#'
#' Basically, only three bits of info are needed from every electrofishing 
#' event: The total catch, T, and a kind of cumulative catch measure termed X
#' in Seber (1978)
#'
#' The data.frame must have the number of passes stored in a colum called 'Runs'
#' And fishing area must be in a row called 'Area'
#'
#' @param data a data.frame containing all relavent info.  
#' @param passnames a character vector giving the names of the columns in which
#'                  the catch data reside.  These must be ordered so that the 
#'                  first comes first, etc. 
#' @return a data frame
#' @export
#' @examples
#' # none yet
getData <- function(data, passnames = paste0("S0_R", 1:6)) {
  catch <- as.matrix(data[passnames])
  rownames(catch) <- NULL
  s <- data $ Runs

  T <- sapply(1:nrow(catch), function(i) sum(catch[i,1:s[i]]))
  maxs <- ncol(catch)
  X <- colSums(sapply(s, function(s) c(s - 1:s, rep(0, maxs - s))) * t(catch), na.rm = TRUE)
  Z <- X / T
  phi <- 1 - Z/(s-1)

  data.frame(
      s = s,
      T = T,
      X = X,
      Z = ifelse(T > 0, Z, 0),
      phi = ifelse(T > 0, phi, 0)
    )
}

#' Utility function to convert parameters to probabilities 
#'
#' The matrix G shoudl be of dimension n x p, 
#' and the parameter vector should be lenght p
#'
#' @param par fitted model parameters
#' @param G The design matrix for a model 
#' @return a data frame
#' @export
#' @examples
#' # none yet
transpar <- function(par, G) {
   1/(1 + exp(-c(G %*% par)))
}
