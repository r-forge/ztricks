\name{Ztricks}
\alias{Ztricks}
\title{Fit an lmer or glmer model with extended random effects}
\description{
  The \pkg{lme4} package provides for specification of random-effects
  terms in linear and generalized linear mixed models using a formula
  language, which is not always adequate.  The \code{Ztricks} function
  allow custom specification of components of the random effects in
  terms of the model matrix.
}
\usage{
Ztricks(formula, data, family = NULL, Znew = list(),
        REML = TRUE, start = NULL, nAGQ = 1, verbose = FALSE,
        subset, weights, na.action, offset, contrasts = NULL,
        model = TRUE, x = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Znew}{a list with one component for each term to be added to the
    random effects.  Each component should itself be a list with
    components named \code{fl}, a named list of the grouping factor,
    \code{Zt} the sparse matrix representation (class
    \code{"\linkS4class{dgCMatrix}"}) of the transpose
    of the columns to be added to \code{Z} and \code{ST}, a dense matrix
    to be added to the \code{ST} list.}
  \item{formula, data, family, REML, start, nAGQ,
    verbose, subset, weights, na.action, offset,
    contrasts, model, x, \dots}{Arguments for the corresponding call to
    \code{\link[lme4]{lmer}} or \code{\link[lme4]{glmer}} in the \pkg{lme4} package.}
}
%\details{}
\value{
  An object of S4 class \code{"\linkS4class{mer}"}.
}
%\references{}
%\author{}
%\note{}
\seealso{\code{link[lme4]{lmer}}, \code{link[lme4]{glmer}}}
\examples{
data(sleepstudy, package = "lme4")
mm <- as(sleepstudy$Subject,"sparseMatrix")
mm@x <- as.double(sleepstudy$Days)

Ztricks(Reaction~Days+(1|Subject),sleepstudy,
  Znew=list(list(
     fl=list(Subject=sleepstudy$Subject),
     Zt=mm,
     ST=matrix(0,1,1,dimnames=list("Days","Days"))
     )))

## This is equivalent to the fit
lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
