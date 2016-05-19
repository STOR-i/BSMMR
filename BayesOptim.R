library(DiceOptim, quietly = T)
setwd("./Results")

## Read the already considered omega and associated MSE from the exteral text file
obs   <- read.table("cross_validation.txt")
obs   <- as.matrix(obs)
omega <- obs[,1]
MSE   <- apply(obs[,-1], 1, mean)
Var   <- apply(obs[,-1], 1, var)

if(length(omega)>30) write.table(c(0,0), "proposal.txt", quote=F, row.names=F, col.names=F, sep=" ")

## Fit Gaussian process with uncertainty in observations
model     <- km(design=data.frame(omega), response = data.frame(MSE), noise.var = Var), 
                control=list(pop.size=50,trace=FALSE), covtype="gauss")

## Calculate the expected improvement for the considered interval 
x         <- seq(0, max(omega), length.out = 1000)
EI_values <- apply(as.matrix(x), 1, EI, model, type="SK")

##Derive the omega resulting in the highest expected improvement
est       <- which(EI_values==max(EI_values))

## Write proposed omega and the corresponding expected improvement to external file
write.table(c(x[est], EI_values[est]), "proposal.txt", quote=F, row.names=F, col.names=F, sep=" ")
if(length(omega)>30) write.table(c(0,0), "proposal.txt", quote=F, row.names=F, col.names=F, sep=" ")

