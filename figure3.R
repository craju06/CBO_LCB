# install.packages(c("tgp","laGP"))
library(tgp)
library(laGP)
options(warn=-1)

# The objective function
obj = function(x1,x2){
  value = -(cos((x1-.1)*x2))^2 - x1*sin(3*x1+x2)
  return(value)
}

# The constraint function
con1 = function(x1,x2){
  t = atan2(x1,x2)
  value = x1^2 + x2^2 -((2*cos(t)-1/2*cos(2*t)-1/4*cos(3*t)-1/8*cos(4*t))^2) - ((2*sin(t))^2)
  return(value)
}

# Objective and constraint function for use 
# with laGP package optimization algorithms
modtom <- function (X, known.only = FALSE)
{
  if(is.null(nrow(X))) X <- matrix(X, nrow=1)
  f <- -(cos((X[,1]-.1)*X[,2]))^2 - X[,1]*sin(3*X[,1]+X[,2])
  #if(known.only) stop("objective is not known")
  t <- atan2(X[,1], X[,2])
  c1 <- X[,1]^2 + X[,2]^2 
  c1 <- c1 - (2*cos(t) - 0.5*cos(2*t) - 0.25*cos(3*t) - 0.125*cos(4*t))^2
  c1 <- c1 - (2*sin(t))^2
  return(list(obj=f, c=cbind(c1)))
}

# S = 100 # The number of additional inputs to select
# n.results = 1 # Number of Monte Carlo runs
# init_input = 20 # Number of initial sample

cat("!! All the arguments must be integer.\n")
cat("!! Number of additional sample must be grater than number of initial sample (>=10).\n")

# Read input from stdin
read_input <- function(default_input){
  input <- file("stdin", "r")
  input <- scan(file=input, what=character(), n=1, quiet=TRUE)
  if(input == 0){
    return(default_input)
  }else{
    return(as.integer(input))
  }
}

cat("Enter Number of initial sample (input 0 for default 20): ")
init_input <- read_input(20)
cat("Enter Number of additional sample (input 0 for default 100): ")
S <- read_input(100)
cat("Enter Number of Monte Carlo runs (input 0 for default 30): ")
n.results <- read_input(30)
############################################################################
### Log-Barrier method using 1/sigma^2 for gamma
############################################################################
print("Log-Barrier method, OOSS:")
RESULTS.f = matrix(NA,nrow=n.results,ncol=(init_input+S))
RESULTS.c1 = matrix(NA,nrow=n.results,ncol=(init_input+S))
RESULTS.x = matrix(NA,nrow=n.results,ncol=(init_input+S))

for(result.iter in 1:n.results){
  print(paste("Simulation:",result.iter))
  
  set.seed(result.iter) # Set a seed for reproducibility purpose
  
  x = lhs(n=init_input,rect=rbind(c(-2.25,2.5),c(-2.5,1.75))) # Initial Latin hypercube sample
  z = obj(x[,1],x[,2])
  w = con1(x[,1],x[,2])
  xx = lhs(n=1000,rect=rbind(c(-2.25,2.5),c(-2.5,1.75)))
  
  for(i in 1:S){
    
    model.f = laGP(xx, 6, init_input-2+i, x, z)   # Surrogate for the objective function
    model.c1 = laGP(xx, 6, init_input-2+i, x, w)  # Surrogate for the constraint function
    
    mu.f =  model.f$mean
    mu.c1 = model.c1$mean
    
    sigma.f = model.f$s2
    sigma.c1 = model.c1$s2
    
    gamma = 1/sigma.f
    E = mu.f - (1/gamma)*((log(-mu.c1) + sigma.c1/(2*mu.c1^2)))  # The acquisition function
    idx = which.min(E)
    
    x = rbind(x,xx[idx,])
    z = c(z,obj(xx[idx,1],xx[idx,2]))
    w = c(w,con1(xx[idx,1],xx[idx,2]))
    
    xx = lhs(n=1000,rect=rbind(c(-2.25,2.5),c(-2.5,1.75)))
    # print(paste(i,result.iter))
    
  }
  
  RESULTS.f[result.iter,] = z
  RESULTS.c1[result.iter,] = w
  
}

# Calculate the best feasible objective value over the function evalutions 
fx = RESULTS.f
for(i in 1:nrow(fx)){
  
  idx = which(RESULTS.c1[i,]<=0)
  if(idx[1] != 1){
    fx[i,1:(idx[1]-1)] = Inf
  }
  hx = rep(0,ncol(fx))
  hx[idx] = RESULTS.f[i,idx]
  
  for(j in 1:(ncol(fx)-1)){
    
    if(RESULTS.f[i,j+1]<fx[i,j] & hx[j+1]!=0 ){
      fx[i,j+1] = RESULTS.f[i,j+1]
    }else{
      fx[i,j+1] = fx[i,j]
    }
    
  }
  
}


# Calculate the average and quantiles
E.f5 = apply(fx,2,mean)
lower5 = apply(fx,2,quantile,prob=0.05)
upper5 = apply(fx,2,quantile,prob=0.95)

OOSS_F <- E.f5

# c(E.f5[(init_input+S)*0.4],E.f5[(init_input+S)*0.65],E.f5[init_input+S])
# c(lower5[(init_input+S)*0.4],lower5[(init_input+S)*0.65],lower5[init_input+S])
# c(upper5[(init_input+S)*0.4],upper5[(init_input+S)*0.65],upper5[init_input+S])

# min_OOSS = matrix(NA,nrow=nrow(RESULTS.f))
# for(i in 1:nrow(RESULTS.f)){
#   idx = which(RESULTS.c1[i,]<=0)
#   min_OOSS[i] = min(RESULTS.f[i,][idx])
# }

infeas_ooss <- matrix(NA,nrow=init_input+S)
for (j in 1:(init_input+S)){
  count = 0
  for (i in 1:n.results){
    if(RESULTS.c1[i,j] >0){
      count = count+1
    }
  }
  infeas_ooss[j] = count
}

############################################################################
### Log-Barrier method using EI-OOSS in place of f(x)
############################################################################
print("Log-Barrier method, EI-OOSS:")
RESULTS.f = matrix(NA,nrow=n.results,ncol=(init_input+S))
RESULTS.c1 = matrix(NA,nrow=n.results,ncol=(init_input+S))
RESULTS.x = matrix(NA,nrow=n.results,ncol=(init_input+S))

for(result.iter in 1:n.results){
  print(paste("Simulation:",result.iter))
  set.seed(result.iter)
  
  x = lhs(n=init_input,rect=rbind(c(-2.25,2.5),c(-2.5,1.75)))
  z = obj(x[,1],x[,2])
  w = con1(x[,1],x[,2])
  xx = lhs(n=1000,rect=rbind(c(-2.25,2.5),c(-2.5,1.75)))
  
  for(i in 1:S){
    
    model.f = laGP(xx, 6, init_input-2+i, x, z)
    model.c1 = laGP(xx, 6, init_input-2+i, x, w)
    
    mu.f =  model.f$mean
    mu.c1 = model.c1$mean
    
    sigma.f = model.f$s2
    sigma.c1 = model.c1$s2
    
    gamma = 1/sigma.f
    fmin = min(z[w<=0])
    E = -((fmin-mu.f)*pnorm((fmin-mu.f)/sqrt(sigma.f)) + sqrt(sigma.f)*dnorm((fmin-mu.f)/sqrt(sigma.f))) - (1/gamma)*((log(-mu.c1) + sigma.c1/(2*mu.c1^2)))
    
    idx = which.min(E)
    
    x = rbind(x,xx[idx,])
    z = c(z,obj(xx[idx,1],xx[idx,2]))
    w = c(w,con1(xx[idx,1],xx[idx,2]))
    
    xx = lhs(n=1000,rect=rbind(c(-2.25,2.5),c(-2.5,1.75)))
    # print(paste(i,result.iter))
    
  }
  
  RESULTS.f[result.iter,] = z
  RESULTS.c1[result.iter,] = w
  
}

# Calculate the best feasible objective value over the function evalutions 
fx = RESULTS.f
for(i in 1:nrow(fx)){
  
  idx = which(RESULTS.c1[i,]<=0)
  if(idx[1] != 1){
    fx[i,1:(idx[1]-1)] = Inf
  }
  hx = rep(0,ncol(fx))
  hx[idx] = RESULTS.f[i,idx]
  
  for(j in 1:(ncol(fx)-1)){
    
    if(RESULTS.f[i,j+1]<fx[i,j] & hx[j+1]!=0 ){
      fx[i,j+1] = RESULTS.f[i,j+1]
    }else{
      fx[i,j+1] = fx[i,j]
    }
    
  }
  
}


# Calculate the average and quantiles
E.f6 = apply(fx,2,mean)
lower6 = apply(fx,2,quantile,prob=0.05)
upper6 = apply(fx,2,quantile,prob=0.95)

EIOOSS_F <- E.f6

# c(E.f6[(init_input+S)*0.4],E.f6[(init_input+S)*0.65],E.f6[init_input+S])
# c(lower6[(init_input+S)*0.4],lower6[(init_input+S)*0.65],lower6[init_input+S])
# c(upper6[(init_input+S)*0.4],upper6[(init_input+S)*0.65],upper6[init_input+S])


# min_EIOOSS = matrix(NA,nrow=nrow(RESULTS.f))
# for(i in 1:nrow(RESULTS.f)){
#   idx = which(RESULTS.c1[i,]<=0)
#   min_EIOOSS[i] = min(RESULTS.f[i,][idx])
# }

infeas_eiooss <- matrix(NA,nrow=init_input+S)
for (j in 1:(init_input+S)){
  count = 0
  for (i in 1:n.results){
    if(RESULTS.c1[i,j] >0){
      count = count+1
    }
  }
  infeas_eiooss[j] = count
}

############################################################################
### Augmented Lagrangian (AL) approach
############################################################################
print("Augmented Lagrangian method, AL:")
out.AL = matrix(NA,nrow=n.results,ncol=(init_input+S))
B = matrix(c(-2.25, -2.5, 2.5, 1.75), ncol=2)
out.AL.C = matrix(NA,nrow=n.results,ncol=(init_input+S))

for(result.iter in 1:n.results){
  print(paste("Simulation:",result.iter))
  set.seed(result.iter)
  x = lhs(n=init_input,rect=rbind(c(-2.25,2.5),c(-2.5,1.75)))
  out.AL.F = optim.auglag(modtom, B, fhat = TRUE, end=init_input+S, Xstart = x,verb = 0)
  out.AL[result.iter,] = out.AL.F$prog
  out.AL.C[result.iter,] = out.AL.F$C
}

E.f7 = apply(out.AL,2,mean)
lower7 = apply(out.AL,2,quantile,prob=0.05)
upper7 = apply(out.AL,2,quantile,prob=0.95)

AL_F <- E.f7

# c(E.f7[(init_input+S)*0.4],E.f7[(init_input+S)*0.65],E.f7[(init_input+S)])
# c(lower7[(init_input+S)*0.4],lower7[(init_input+S)*0.65],lower7[(init_input+S)])
# c(upper7[(init_input+S)*0.4],upper7[(init_input+S)*0.65],upper7[(init_input+S)])

# min_AL = matrix(NA,nrow=n.results)
# for (i in 1:n.results){
#   min_AL[i] = min(out.AL[i,])
# }

infeas_al <- matrix(NA,nrow=(init_input+S))
for (j in 1:(init_input+S)){
  count = 0
  for (i in 1:n.results){
    if(out.AL.C[i,j] >0){
      count = count+1
    }
  }
  infeas_al[j] = count
}

############################################################################
### Constrained expected improvement (CEI) approach
############################################################################
print("Constrained expected improvement method, CEI:")
out.CEI = matrix(NA,nrow=n.results,ncol=(init_input+S))
B = matrix(c(-2.25, -2.5, 2.5, 1.75), ncol=2)
out.CEI.C = matrix(NA,nrow=n.results,ncol=(init_input+S))

for(result.iter in 1:n.results){
  print(paste("Simulation:",result.iter))
  set.seed(result.iter)
  x = lhs(n=init_input,rect=rbind(c(-2.25,2.5),c(-2.5,1.75)))
  out.CEI.F = optim.efi(modtom, B, fhat = TRUE, end=init_input+S, Xstart = x,verb = 0)
  out.CEI[result.iter,] = out.CEI.F$prog
  out.CEI.C[result.iter,] = out.CEI.F$C
}

E.f8 = apply(out.CEI,2,mean)
lower8 = apply(out.CEI,2,quantile,prob=0.05)
upper8 = apply(out.CEI,2,quantile,prob=0.95)

CEI_F <- E.f8

# c(E.f8[(init_input+S)*0.4],E.f8[(init_input+S)*0.65],E.f8[(init_input+S)])
# c(lower8[(init_input+S)*0.4],lower8[(init_input+S)*0.65],lower8[(init_input+S)])
# c(upper8[(init_input+S)*0.4],upper8[(init_input+S)*0.65],upper8[(init_input+S)])

# min_CEI = matrix(NA,nrow=n.results)
# for (i in 1:n.results){
#   min_CEI[i] = min(out.CEI[i,])
# }

infeas_cei <- matrix(NA,nrow=(init_input+S))
for (j in 1:(init_input+S)){
  count = 0
  for (i in 1:n.results){
    if(out.CEI.C[i,j] >0){
      count = count+1
    }
  }
  infeas_cei[j] = count
}
############################################################################
### End of algorithms
############################################################################

## comment out
command1 <- paste("python3 figure3.py", init_input, S, n.results)
command2 <- paste("python figure3.py", init_input, S, n.results)
is_python3 <- system("python3 --version", intern = TRUE) != ""
if (length(is_python3)!=0) {
  print("Its python3")
  system(command1)
  } else {
  print("Its python")
  system(command2)
}
## comment out


MLCB <- read.csv('MLCB.csv',header=F) # M-based lower confidence bound
MLCB_F <- MLCB$V1
file.remove("MLCB.csv")

infeas_mlcb <- read.csv('infeas_mlcb.csv',header=F) # M-based lower confidence bound
infeas_mlcb <- infeas_mlcb$V1
file.remove("infeas_mlcb.csv")
#######################################################
####### Final Plot
####################################################

jpeg("Plot.jpeg", width = 12, height = 12, units = "cm", res = 600)
cex = 0.7
par(mfrow = c(2, 1),mar = c(4, 4, 1, 1),cex.lab=cex, cex.axis=cex, cex.main=cex)
# par(mfrow = c(2, 1),mar = c(4, 4, .5, .5))
plot(MLCB_F,type="l",lwd=2,xlab="Blackbox evaluations (n)", 
     ylab="Best valid objective (f)",axes=FALSE,ylim=c(-2.05,-1), lty = 1,col=1)
lines(OOSS_F,col=4,lwd=2, lty = 4)
lines(EIOOSS_F,col=5,lwd=2, lty = 5)
lines(AL_F,col=6,lwd=2, lty = 6)
lines(CEI_F,col=7,lwd=2, lty = 7)
axis(1)
axis(2)
legend("topright",c("M-LCB","OOSS","EI-OOSS","AL","CEI","Global Solution"),
       lty=c(1,4,5,6,7,2),lwd=c(2,2,2,2,2,1),bty="n",
       col=c(1,4,5,6,7,2),cex=cex)
abline(h=obj(2.0052938,1.1944509),lty=2,lwd=1,col=2)

plot(infeas_mlcb,type="l",lwd=2,xlab="Blackbox evaluations (n)",
     ylab="Number of infeasible points",axes=FALSE,ylim=c(0,n.results),lty = 1,col=1)
lines(infeas_ooss,col=4,lwd=2, lty = 4)
lines(infeas_eiooss,col=5,lwd=2, lty = 5)
lines(infeas_al,col=6,lwd=2, lty = 6)
lines(infeas_cei,col=7,lwd=2, lty = 7)
axis(1)
axis(2)
legend("topleft",c("M-LCB","OOSS","EI-OOSS","AL","CEI","Zero line"),
       lty=c(1,4,5,6,7,2),lwd=c(2,2,2,2,2,1),bty="n",
       col=c(1,4,5,6,7,2),cex=cex)
abline(h=0.0,lty=2,lwd=1,col=2)
dev.off()

system("open Plot.jpeg")