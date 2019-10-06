library(ggplot2) 
library(dplyr) 
library(gridExtra) 
set.seed(1)
library(MASS)

#####Part C#####
# generate data
n          =     100
p          =     2
X          =     matrix(ncol = p, nrow = n)

p1         =     0.6
p2        =    .3
p3        =     .1

rho1         =     0.25
rho2        =    .5
rho3        =     .75
o         =     0.5

C1         =     o^2 * matrix(c(1, rho1, rho1, 1), nrow = 2)
C2         =     o^2 * matrix(c(1, rho2, rho2, 1), nrow = 2)
C3         =     o^2 * matrix(c(1, rho3, rho3, 1), nrow = 2)

mu1        =     c(0,1)
mu2        =     c(1,0)
mu3        =     c(-1,1)


# Data generating process
y          =     as.factor(sample(3, n, replace = TRUE, prob = c(p1,p2,p3)))
for (i in 1:n){
  mu     =     (y[i]==1)*mu1 + (y[i]==2)*mu2 + (y[i]==3)*mu3
  C      =     (y[i]==1)*C1 + (y[i]==2)*C2 + (y[i]==3)*C3
  X[i,]  =     mvrnorm(1, mu, C)
}

dataf     =     data.frame(y, X)
plot1     =     ggplot(dataf, aes(x=X1, y=X2, colour=y))+geom_point()+ggtitle("scatter plot")
grid.arrange(plot1, nrow =1)

########Part D############
# separate the features corresponding to 1, 2, and 3
n1        =   sum(y==1)
n2        =   sum(y==2)
n3        =   sum(y==3)
X1        =   X[y==1,]
X2        =   X[y==2,]
X3        =   X[y==3,]
mu1.hat   =   colMeans(X1)
mu2.hat   =   colMeans(X2)
mu3.hat   =   colMeans(X3)
p1.hat    =   n1/n
p2.hat    =   n2/n
p3.hat    =   n3/n


X.grid     =     data.matrix(expand.grid(x=seq(from = min(X[,1]), to =  max(X[,1]), length.out = 100),
                                         y=seq(from = min(X[,2]), to =  max(X[,2]), length.out = 100)))

qda_discrim = function(x_, mu_, cov_, p_) {
  d = -0.5 * t(x_ - mu_) %*% solve(cov_) %*% (x_ - mu_) - 0.5 * log( det(cov_) )  + log(p_)
  return(d)
}

#results of discriminant function stored as d1, d2, d3 where di = results of function

d1 = apply(X.grid, 1, qda_discrim, mu_ = mu1, cov_ = C1, p_=.6)
d2 = apply(X.grid, 1, qda_discrim, mu_ = mu2, cov_ = C2, p_=.3)
d3 = apply(X.grid, 1, qda_discrim, mu_ = mu3, cov_ = C3, p_=.1)

y.hat.grid = ifelse((d1>d2) & (d1>d3), 1, ifelse((d2>d1) & (d2>d3), 2, 3))

dataff     =     data.frame(X.grid, y.hat.grid)
dataf      =     data.frame(y, X)

p1         =     ggplot(dataf, aes(x=X1, y=X2, colour=y))+geom_point()+ggtitle("Data")
p2         =     ggplot(dataff, aes(x=x, y=y, colour=as.factor(y.hat.grid)))+geom_point()+scale_color_discrete()+ggtitle("Optimal QDA")
grid.arrange(p1,  p2, nrow =1) 

#### Part E ##############

d1 = apply(X.grid, 1, qda_discrim, mu_ = mu1.hat, cov_ = cov(X1), p_=p1.hat)
d2 = apply(X.grid, 1, qda_discrim, mu_ = mu2.hat, cov_ = cov(X2), p_=p2.hat)
d3 = apply(X.grid, 1, qda_discrim, mu_ = mu3.hat, cov_ = cov(X3), p_=p3.hat)

y.hat.grid = ifelse((d1>d2) & (d1>d3), 1, ifelse((d2>d1) & (d2>d3), 2, 3))

dataff     =     data.frame(X.grid, y.hat.grid)

ggplot(dataff, aes(x=x, y=y, colour=as.factor(y.hat.grid)))+geom_point()+scale_color_discrete()+ggtitle("QDA")

#### Part F1 ###############
# generate data
n          =     1000
p          =     2
X          =     matrix(ncol = p, nrow = n)

p1         =     0.6
p2        =    .3
p3        =     .1

rho1         =     0.25
rho2        =    .5
rho3        =     .75
o         =     0.5

C1         =     o^2 * matrix(c(1, rho1, rho1, 1), nrow = 2)
C2         =     o^2 * matrix(c(1, rho2, rho2, 1), nrow = 2)
C3         =     o^2 * matrix(c(1, rho3, rho3, 1), nrow = 2)

mu1        =     c(0,1)
mu2        =     c(1,0)
mu3        =     c(-1,1)


# Data generating process
y          =     as.factor(sample(3, n, replace = TRUE, prob = c(p1,p2,p3)))
for (i in 1:n){
  mu     =     (y[i]==1)*mu1 + (y[i]==2)*mu2 + (y[i]==3)*mu3
  C      =     (y[i]==1)*C1 + (y[i]==2)*C2 + (y[i]==3)*C3
  X[i,]  =     mvrnorm(1, mu, C)
}

dataf     =     data.frame(y, X)
plot1     =     ggplot(dataf, aes(x=X1, y=X2, colour=y))+geom_point()+ggtitle("scatter plot")
grid.arrange(plot1, nrow =1)
n1        =   sum(y==1)
n2        =   sum(y==2)
n3        =   sum(y==3)
X1        =   X[y==1,]
X2        =   X[y==2,]
X3        =   X[y==3,]
mu1.hat   =   colMeans(X1)
mu2.hat   =   colMeans(X2)
mu3.hat   =   colMeans(X3)
p1.hat    =   n1/n
p2.hat    =   n2/n
p3.hat    =   n3/n


X.grid     =     data.matrix(expand.grid(x=seq(from = min(X[,1]), to =  max(X[,1]), length.out = 100),
                                         y=seq(from = min(X[,2]), to =  max(X[,2]), length.out = 100)))

#results of discriminant function stored as d1, d2, d3 where di = results of function

#### Part F2 ##############

d1 = apply(X.grid, 1, qda_discrim, mu_ = mu1, cov_ = C1, p_=.6)
d2 = apply(X.grid, 1, qda_discrim, mu_ = mu2, cov_ = C2, p_=.3)
d3 = apply(X.grid, 1, qda_discrim, mu_ = mu3, cov_ = C3, p_=.1)

y.hat.grid = ifelse((d1>d2) & (d1>d3), 1, ifelse((d2>d1) & (d2>d3), 2, 3))

dataff     =     data.frame(X.grid, y.hat.grid)
dataf      =     data.frame(y, X)

p1         =     ggplot(dataf, aes(x=X1, y=X2, colour=y))+geom_point()+ggtitle("Data")
p2         =     ggplot(dataff, aes(x=x, y=y, colour=as.factor(y.hat.grid)))+geom_point()+scale_color_discrete()+ggtitle("Optimal QDA")
grid.arrange(p1,  p2, nrow =1) 

#### Part F3 ##############

d1 = apply(X.grid, 1, qda_discrim, mu_ = mu1.hat, cov_ = cov(X1), p_=p1.hat)
d2 = apply(X.grid, 1, qda_discrim, mu_ = mu2.hat, cov_ = cov(X2), p_=p2.hat)
d3 = apply(X.grid, 1, qda_discrim, mu_ = mu3.hat, cov_ = cov(X3), p_=p3.hat)

y.hat.grid = ifelse((d1>d2) & (d1>d3), 1, ifelse((d2>d1) & (d2>d3), 2, 3))

dataff     =     data.frame(X.grid, y.hat.grid)

ggplot(dataff, aes(x=x, y=y, colour=as.factor(y.hat.grid)))+geom_point()+scale_color_discrete()+ggtitle("QDA")
