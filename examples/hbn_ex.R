library(lineprof)

##############################################
# An implementation of Hub Binary Network
##############################################
set.seed(1000)
n=250
p=30

# A network with 2 hubs
network<-HubNetwork(p,0.95,2,0.1,type="binary")
Theta <- network$Theta
truehub <- network$hubcol
# The four hub nodes have indices 4,5
print(truehub)

# Generate data matrix x
X <- binaryMCMC(n,Theta,burnin=500,skip=100)

# Run Hub Binary Network to estimate Theta
Rprof("L_C.out")
res1 <- hbn(X,2,1,3,trace=TRUE)
Rprof(NULL)
S_C <- summaryRprof("L_C.out")

Rprof("L_R.out")
res1 <- hbn(X,2,1,3,trace=TRUE, method="R")
Rprof(NULL)
S_R <- summaryRprof("L_R.out")
speedup <- S_R$sampling.time / S_C$sampling.time 
speedup
# print out a summary of the object hbn
summary(res1)

# We see that the estimated hub nodes have indices 1,5
# We successfully recover the hub nodes

# Plot the resulting network
#plot(res1) 