##################################################
## simple lift function.
##################################################
#--------------------------------------------------
#plots lift function
liftf = function(yl,phatl,dopl=TRUE) {
if(is.factor(yl)) yl = as.numeric(yl)-1
oo = order(-phatl)
sy = cumsum(yl[oo])/sum(yl==1)
if(dopl) {
   par(mai=c(.8,1.2,.5,.5))
   ii = (1:length(sy))/length(sy)
   plot(ii,sy,type='l',lwd=2,col='blue',xlab='% tried',ylab='% of successes',cex.lab=2)
   abline(0,1,lty=2)
}
return(sy)
}
##################################################
## simple lift function for list of phats.
##################################################
#--------------------------------------------------
#plots lift function
liftfL = function(yl,phatlL) {
if(is.factor(yl)) yl = as.numeric(yl)-1
n=length(yl)
par(mai=c(.8,1.2,.5,.5))
plot(c(0,1),c(0,1),type="n",xlab='% tried',ylab='% of successes',cex.lab=2)
abline(0,1,lty=2)
ii = (1:n)/n
for(i in 1:length(phatlL)) {
   oo = order(-phatlL[[i]])
   sy = cumsum(yl[oo])/sum(yl==1)
   lines(ii,sy,type="l",lwd=2,col=i)
}
}
############################################################
dummies<-function(cl)
{
  n <- length(cl)
  cl <- as.factor(cl)
  x <- matrix(0, n, length(levels(cl)) )
  x[(1:n) + n*(unclass(cl)-1)] <- 1
  dimnames(x) <- list(names(cl), levels(cl))
  as.data.frame(x)
}

CheckBinFit<-function(y,phat,nq=20)
{
if(is.factor(y))
{ 
y<-as.double(y)
}
y<-y-mean(y)
y[y>0]<-1
y[y<=0]<-0
quants<-quantile(phat,probs=(1:nq)/(nq+1))
names(quants)<-NULL
quants<-c(0,quants,1)
phatD<-rep(0,nq+1)
phatF<-rep(0,nq+1)
nv <- rep(0,nq+1)
for(i in 1:(nq+1))
{
which<-((phat<=quants[i+1])&(phat>quants[i]))
nn <- length(phat[which])
if(nn>0){
nv[i] <- nn
phatF[i]<-mean(phat[which])
phatD[i]<-mean(y[which])
} else {
phatF[i] <- -1
phatD[i] <- -1
}
}
plot(phatF[phatF!=-1],phatD[phatD!=-1],xlab="phat",ylab="data")
abline(0,1)
list(phat=phatF,data=phatD,nv=nv)
}


lift<-function(y,pl)
{
n<-length(y)
np<-length(pl)
yy<-as.numeric(y)
yy<-yy-mean(yy)
yyy<-rep(0,n)
yyy[yy>0]<-1
nsuc<-sum(yyy)
liftl<-list()
for(i in 1:np)
{
#print(i)
oo<-order(-pl[[i]])
ll<-cumsum(yyy[oo])
liftl[[i]]<-ll/nsuc
}
prop<-mean(yyy)
#win.graph()
#dev.set(dev.cur())
plot(c(0,1),c(0,1),xlab="% obs",ylab="% success",type="l",col=1,lty=2,lwd=2)
lines(c(0,prop),c(0,1),col=1,lty=2,lwd=2)
lines(c(prop,1),c(1,1),col=1,lty=2,lwd=2)
for(i in 1:np)
{
lines((1:n)/n,liftl[[i]],col=i,lwd=2)
}
list(lift=liftl,y=yyy)
}

#Functions written or modified, R. Adam Molnar.
liftn<-
function (resp,pred,doplot=TRUE,firstplot=TRUE,col=1,high=T)  
{ 
#liftn
#Function written for GSB 41201, U. of Chicago. 
#Computes lift function for binary response resp, based on 
#predictor pred.  Accounts for ties correctly. 
#Only resets plot if firstplot=T, so can be used for multiples. 
#col is the color number for the lift curve plot. 
#high is T if high values of the prediction are to be considered
#successes, and thus picked first.  Set it to F if for some
#reason low values of pred should be ordered first.

#Use sum and cumsum to get cumulative proportions going down the table. 
#This relies on the inherent sort of pred by table(). 
if (mode(pred)=="list") {pred <- pred[[1]]} 
#Sort based on rownames, which are the possible values of pred.
ttab <- table(pred,resp) 
if(high) {ttab <- ttab[rev(1:nrow(ttab)),]}
sum0 <- sum(ttab[,1]) 
sum1 <- sum(ttab[,2]) 
#The leading zero is needed to get nice plots. 
cs0 <- c(0, cumsum(ttab[,1])) / sum0 
cs1 <- c(0, cumsum(ttab[,2])) / sum1 
cst <- c(0, cumsum(ttab[,1]+ttab[,2])) / (sum0+sum1) 
 
#Plotting, if doplot=T. 
if(doplot) 
        { 
        if(firstplot) 
                { 
                plot(c(0,1), c(0,1), type="l", col="gray", lty=2, lwd=1, 
                xlab="% of observations", ylab="% of successes")  
                prop <- sum1 / (sum0+sum1) 
                lines(c(0,prop), c(0,1), col="gray", lty=2, lwd=1) 
                lines(c(prop,1), c(1,1), col="gray", lty=2, lwd=1) 
                } 
        lines(cst, cs1, lwd=2, lty=1, col=col) 
        } 
#SAS C computation is really percentage of comparisons won. 
#This code uses "high" implicitly in the ordering. 
s1 <- c(0,ttab[,2]) 
lt <- nrow(ttab) + 1 
win <- (1 - cs0) + (1 - c(0, cs0[-lt])) 
sasc <- sum(s1 * win) / (2 * sum1)  
return(sasc) 
}


liftl <- 
function(yl,pl,doplot=TRUE)
{
#liftl, function written for GSB 41201, U. of Chicago.
#In keeping with lift() function, takes a list of response vectors as yl,
#and a list of prediction vectors as pl.  Vectors may be of different lengths.
#Uses calls to lift() to process elements.
np <- length(yl)
ccomp <- rep(0,np)
for(i in 1:np)
	{
	if (i==1) {fp <- TRUE} else {fp <- FALSE}
	ccomp[i] <- liftn(yl[[i]], pl[[i]], doplot, fp, i)
	}
return(ccomp)
}

loss <- 
function(y, pred, eps=.0001)
{
#Loss is a function to compute the appropriate loss,
#based on binary, categorical, or numeric y vector,
#a prediction vector pred, and a minimum probability tolerance eps.
#GSB 41201, U. of Chicago.

val <- 0
#Factors first.
if(is.factor(y))
{
  #Trim away a tree result, using the second or higher value column.
  if((length(levels(y))==2) & is.matrix(pred)) pred <- pred[,2]

  if(!is.matrix(pred))
	{
 	#assume y is 0 or 1 and pred is prob of a 1
	ny <- as.numeric(y)-1
 	pred[pred < eps] <- eps
	pred[pred > (1-eps)] <- 1-eps
	temp <- -2*(ny*log(pred) + (1-ny)*log(1-pred))
 	val <- sum(temp)
	}

  else #Categorical y, loop is quite slow.
	{
	for(i in 1:length(y))
	  {
	  temp <- pred[i,as.character(y[i])]
	  ltemp <- -2*log(max(temp,eps))
	  val<- val+ltemp
	  }
	}
  }

else #Numeric Y
  {val <- sum((y-pred)^2)}

return(val)
}

#doTrainVal = function(xtrain,ytrain,xval,yval,predf,lossf,parm,parc=NULL)
##xtrain: dataframe of train x
##ytrain: vector of train y
##xval: dataframe of val x
##yval: vector of val y
##predf: takes (xtrain,ytrain,xval,parc,parm) and a row of parmspred to make prediction
##lossf: gets loss between output of predf and yval
##parmsfit: parameters needed for fit on (xtrain,ytrain) 
##parc: constant set of parameters used in predicting
##parm: matrix of prediction parameter, ith row is ith setting for evaluation
#{
#}

tree2s<-function(df1,df2,yi,xi,sizev,mindv=1e-6,dpl=TRUE)
{
#setup up in=1 and out =2 data frames
dodf1<-data.frame(y=df1[[yi]])
for(i in 1:length(xi)) dodf1[[i+1]]<-df1[[xi[i]]]
names(dodf1)[2:ncol(dodf1)]<-names(df1)[xi]
dodf2<-data.frame(y=df2[[yi]])
for(i in 1:length(xi)) dodf2[[i+1]]<-df2[[xi[i]]]
names(dodf2)[2:ncol(dodf2)]<-names(df2)[xi]
# get the losses
bigtree<-tree(y~.,data=dodf1,mindev=mindv)
maxsz<-length(unique(bigtree$where))
sizev<-sizev[sizev<=maxsz]
sizev<-sizev[sizev>1]
nt<-length(sizev)
bigtree<-prune.tree(bigtree,best=max(sizev))
inlossv<-rep(0,nt)
outlossv<-rep(0,nt)
for(i in 1:nt){
temp<-prune.tree(bigtree,best=sizev[i])
fitin<-predict(temp,dodf1)
fitout<-predict(temp,dodf2)
inlossv[i]<-loss(dodf1$y,fitin)
outlossv[i]<-loss(dodf2$y,fitout)
}
inlossv<-sqrt(inlossv/nrow(dodf1))
outlossv<-sqrt(outlossv/nrow(dodf2))
if(dpl){
plot(range(sizev),range(c(inlossv,outlossv)),xlab="tree size",ylab="loss",type="n")
points(sizev,inlossv,col="red",pch="i")
points(sizev,outlossv,col="blue",pch="o")
}
list(inloss=inlossv,outloss=outlossv,size=sizev)
}

treeCV<-function(df,yi,xi,tsize,mindv=1e-6,np=10,dpl=TRUE)
{
#first set up group indication ****************************
nob<-nrow(df)
cvsize <- trunc(nob/np)
CvInd <- rep(1:np,rep(cvsize,np))
if(cvsize*np!=nob) CvInd[(np*cvsize+1):nob] <- np
#copy data frame and permute it **************************
dodf<-data.frame(y=df[[yi]])
for(i in 1:length(xi)) dodf[[i+1]]<-df[[xi[i]]]
names(dodf)[2:ncol(dodf)]<-names(df)[xi]
dodf<-dodf[sample(1:nob,nob),]
# loop over quantiles and accumulate results ***********
nts<-length(tsize)
lossV<-rep(0,nts)
for(i in 1:np){
#print(i)
indf<-dodf[CvInd!=i,]
outdf<-dodf[CvInd==i,]
bigTree<-tree(y~.,data=indf,control=tree.control(nobs=nrow(indf),minsize=2,mindev=mindv))
bigTree<-prune.tree(bigTree,best=max(tsize))
for(j in 1:nts){
#print(j)
temptree<-prune.tree(bigTree,best=tsize[j])
temppred<-predict(temptree,outdf)
lossV[j]<-lossV[j] + loss(outdf$y,temppred)
}
}
if(is.factor(dodf$y)){
lossV<-lossV/nob
}
else {
lossV<-sqrt(lossV/nob)
}
if(dpl){
plot(tsize,lossV,xlab="tree size",ylab="cv loss")
}
data.frame(tsize=tsize,loss=lossV)
}
############################################################
############################################################
## Function to do cross validation.
## docv is a general method that takes a prediction function as an argument.
## docvknn is for kNN, it calls docv handing in a wrapper of kknn.
############################################################
############################################################
#--------------------------------------------------
if(1) {#set fit and loss functions
mse=function(y,yhat) {return(sum((y-yhat)^2))}
doknn=function(x,y,xp,k) {
   kdo=k[1]
   train = data.frame(x,y=y)
   test = data.frame(xp); names(test) = names(train)[1:(ncol(train)-1)]
   near  = kknn(y~.,train,test,k=kdo,kernel='rectangular')
   return(near$fitted)
}
}
#--------------------------------------------------
if(1) {#make docv
docv = function(x,y,set,predfun,loss,nfold=10,doran=TRUE,verbose=TRUE,...)
{
   #a little error checking
   if(!(is.matrix(x) | is.data.frame(x))) {cat('error in docv: x is not a matrix or data frame\n'); return(0)}
   if(!(is.vector(y))) {cat('error in docv: y is not a vector\n'); return(0)}
   if(!(length(y)==nrow(x))) {cat('error in docv: length(y) != nrow(x)\n'); return(0)}

   nset = nrow(set); n=length(y) #get dimensions
   if(n==nfold) doran=FALSE #no need to shuffle if you are doing them all.
   cat('in docv: nset,n,nfold: ',nset,n,nfold,'\n')
   lossv = rep(0,nset) #return values
   if(doran) {ii = sample(1:n,n); y=y[ii]; x=x[ii,,drop=FALSE]} #shuffle rows

   fs = round(n/nfold) # fold size
   for(i in 1:nfold) { #fold loop
      bot=(i-1)*fs+1; top=ifelse(i==nfold,n,i*fs); ii =bot:top
      if(verbose) cat('on fold: ',i,', range: ',bot,':',top,'\n')
      xin = x[-ii,,drop=FALSE]; yin=y[-ii]; xout=x[ii,,drop=FALSE]; yout=y[ii]
      for(k in 1:nset) { #setting loop
         yhat = predfun(xin,yin,xout,set[k,],...)
         lossv[k]=lossv[k]+loss(yout,yhat)
      } 
   } 

   return(lossv)
}
}
#cv version for knn
docvknn = function(x,y,k,nfold=10,doran=TRUE,verbose=TRUE) {
return(docv(x,y,matrix(k,ncol=1),doknn,mse,nfold=nfold,doran=doran,verbose=verbose))
}
#--------------------------------------------------
if(0) {#simulate data
set.seed(99)
n=10000
n=500
x=matrix(sort(rnorm(n)),ncol=1) #x must be a matrix (or df?)
y = drop(2*x + .5*rnorm(n)) #must be a column
xp = x[sort(sample(1:n,floor(n/2))),,drop=0]
}

#--------------------------------------------------
if(0) {#test doknn
yhatknn = doknn(x,y,xp,c(50))
plot(x,y)
lines(xp,yhatknn,col='red')
}
#--------------------------------------------------
if(0) {
kv = 2:30
set.seed(99)
tm= system.time({temp1 = docv(x,y,matrix(kv,ncol=1),doknn, mse,nfold=5)})
temp2 = docv(x,y,matrix(kv,ncol=1),doknn, mse,nfold=5)
temp3=docvknn(x,y,kv,nfold=n,verbose=FALSE)
temp4=docvknn(x,y,kv,nfold=10)
ry = range(c(temp1,temp2,temp3))
par(mfrow=c(1,2))
plot(log(1/kv),temp1,type='l',ylim=ry)
lines(log(1/kv),temp2,col='red')
lines(log(1/kv),temp3,col='blue')
lines(log(1/kv),temp4,col='green')
plot(kv,temp1,type='l',ylim=ry)
lines(kv,temp2,col='red')
lines(kv,temp3,col='blue')
lines(kv,temp4,col='green')
}
