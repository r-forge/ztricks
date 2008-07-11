## Check the format of the additional Zt matrix slabs for consistency
chkLen <- function(Znewi, n)
{
    fl <- Znewi$fl
    stopifnot(is.list(fl), length(fl) == 1, length(names(fl)) == 1, is.factor(fl[[1]])) 
    nlev <- length(levels(fli <- fl[[1]][drop = TRUE]))
    stopifnot(length(fli) == n,
              is(Zti <- Znewi$Zt,"sparseMatrix"),
              ncol(Zti) == n)
    nb <- nrow(Zti)
    STi <- as.matrix(Znewi$ST)
    stopifnot(is.numeric(STi), ncol(STi) == nrow(STi),
              length(colnames(STi)) > 0, ncol(STi) * nlev == nb)
    TRUE
}

## Znew should be a list with one component for each term containing
## fl - a named list containing a factor
## Zt - sparse matrix representation of the transpose of the columns to be added to Z
## ST - dense matrix to be added to the ST list

Ztricks <-
    function(formula, data, family = NULL, Znew = list(),
             REML = TRUE, start = NULL, nAGQ = 1, verbose = FALSE,
             subset, weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, ...)
{
    mc <- match.call()
    if (!length(Znew)) { # call lmer directly
        mc$Znew <- NULL
        mc[[1]] <- as.name("lmer")
        return(eval.parent(mc))
    }
    
    ## Call lmer without Znew and with doFit = FALSE
    mc$Znew <- NULL
    mc$doFit <- FALSE
    mc[[1]] <- as.name("lmer")
    lf <- eval.parent(mc)

    ## Check the elements of Znew and add them to the FL list
    n <- length(lf$fr$Y)
    FL <- lf$FL
    ntrm <- length(FL$trms)
    for (i in seq_along(Znew)) {
        Znewi <- Znew[[i]]
        chkLen(Znewi, n)
        FL$trms[[ntrm + i]] <- list(A = Znewi$Zt, Zt = Znewi$Zt, ST = Znewi$ST)
        fl <- Znewi$fl
        
        ## check if the grouping factor has already been used
        ## and adjust the assign attribute accordingly
        if (!(ff <- match(names(fl)[[1]], names(FL$fl), no = 0))) {
            FL$fl <- c(FL$fl, fl)
            attr(lf$FL$fl, "assign") <- c(attr(FL$fl, "assign"), length(FL$fl))
        } else {
            attr(FL$fl, "assign") <- c(attr(FL$fl, "assign"), ff)
        }
    }
    lf$FL <- FL
    ans <- do.call(if (!is.null(lf$glmFit))
                   lme4:::glmer_finalize else lme4:::lmer_finalize, lf)
    ans@call <- match.call()
    ans
}

