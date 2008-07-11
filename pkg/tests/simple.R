library(Ztricks)
data(sleepstudy, package = "lme4")
mm <- as(sleepstudy$Subject,"sparseMatrix")
mm@x <- as.double(sleepstudy$Days)

Ztricks(Reaction~Days+(1|Subject),sleepstudy,
  Znew=list(list(
     fl=list(Subject=sleepstudy$Subject),
     Zt=mm,
     ST=matrix(0,1,1,dimnames=list("Days","Days"))
     )))

