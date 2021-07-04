#################### DM - HW1 ####################

rm(list=ls())
setwd("C:\\Users\\User\\Desktop\\데이터마이닝이론\\data")  
# harris.dat이 존재하는 위치로 지정해주어야한다.



# example

data=readline("Enter the data file name: " )
harris.dat

cat("Select the data coding format(1 = 'a b c' or 2 = 'a,b,c'): ")
fm = scan(n=1, quiet=TRUE)
1
if(fm==1) {form = ""} else {form = ","}
data=read.table(data, sep=form)
data


fit=lm(V1~V2+V3+V4+V5, data=data)
sum = summary(fit)
aov(V1~V2+V3+V4+V5, data=data)



# make multiple linear regression function

regression.model <- function() {
  
  # import data file
  data = readline("Enter the data file name (with extension name) : " )
  cat("Select the data coding format(1 = 'a b c' or 2 = 'a,b,c'): ")
  fm = scan(n=1, quiet=TRUE)
  if(fm==1) {form = ""} else {form = ","}
  data = read.table(data, sep=form)
  
  # design matrix
  n = dim(data)[1]
  p = dim(data)[2]-1
  one = matrix(1, nrow=n, ncol=1)
  I = diag(n)
  
  y = as.matrix(data[,1])
  x = as.matrix(cbind(one,data[,-1]))
  H = x%*%solve(t(x)%*%x)%*%t(x)
  H0 = one%*%solve(t(one)%*%one)%*%t(one)

  # multiple regression result
  b = round(solve(t(x)%*%x)%*%t(x)%*%y, 4)
  yhat = x%*%b
  SST = t(y)%*%(I-H0)%*%y
  SSE = t(y)%*%(I-H)%*%y
  Rsquare = round(1 - SSE/SST, 4)
  MSE = round(SSE/(n-p), 4)
  
  # output - coefficient
  name = paste("Beta",c(0:p),":" ,sep="")
  name[1] = "Constant:"
  row.names(b) = name
  
  # output - ID, Actual values, Fitted values
  y.values = cbind(c(1:n), y, round(yhat,1))
  
  # output- Model Summary
  model.smy = rbind(Rsquare, MSE)
  row.names(model.smy) <- c("R-square = ","MSE = ")
  
  # make output file
  outputname = readline("Write the output file name you want to save (without extension name) : ")
  outputname = paste(outputname,".txt",sep="")

  cat("Coefficients","\n","-------------","\n",file = outputname,sep="")
  write.table(b, outputname, sep= " ", row.names=TRUE, col.names=FALSE, append=TRUE, quote=FALSE)
  
  cat("\n","ID, Actual values, Fitted values","\n","--------------------------------","\n",file = outputname,sep="", append=TRUE)
  write.table(y.values, outputname, sep= ", ", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
  
  cat("\n","Model Summary","\n","-------------","\n",file = outputname,sep="", append=TRUE)
  write.table(model.smy, outputname, sep= "", row.names=TRUE, col.names=FALSE, append=TRUE, quote=FALSE)
  
  cat("Output file has been saved in ",getwd(),"/",outputname,sep="")
}


regression.model()
harris.dat
1
harris.regression.output
