
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
    function(formula, data, Znew = list(), REML = TRUE,
             start = NULL, verbose = FALSE, subset,
             weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, ...)
{
    mc <- match.call()
    if (!length(Znew)) {
        mc$Znew <- NULL
        mc[[1]] <- as.name("lmer")
        return(eval(mc))
    }
    stopifnot(length(formula <- as.formula(formula)) == 3)

    fr <- lme4:::lmerFrames(mc, formula, contrasts) # model frame, X, etc.
    FL <- lme4:::lmerFactorList(formula, fr$mf, 0L, 0L) # flist, Zt, cnames
    Y <- as.double(fr$Y)

    n <- length(Y)
    ntrm <- length(FL$trms)
    browser()
    for (i in seq_along(Znew)) {
        Znewi <- Znew[[i]]
        chkLen(Znewi, n)
        FL$trms[[ntrm + i]] <- list(A = Znewi$Zt, Zt = Znewi$Zt, ST = Znewi$ST)
        fl <- Znewi$fl
        if (!(ff <- match(names(fl)[[1]], names(FL$fl), no = 0))) {
            FL$fl <- c(FL$fl, fl)
            attr(FL$fl, "assign") <- c(attr(FL$fl, "assign"), length(FL$fl))
        } else {
            attr(FL$fl, "assign") <- c(attr(FL$fl, "assign"), ff)
        }
    }
    lme4:::lmer_finalize(mc, fr, FL, start, REML, verbose)
}

mm <- as(sleepstudy$Subject,"sparseMatrix")
mm@x <- as.double(sleepstudy$Days)

Ztricks(Reaction~Days+(1|Subject),sleepstudy,
  Znew=list(list(
     fl=list(Subject=sleepstudy$Subject),
     Zt=mm,
     ST=matrix(0,1,1,dimnames=list("Days","Days"))
     )))


