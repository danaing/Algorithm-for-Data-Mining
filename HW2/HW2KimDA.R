#################### DM - HW2 ####################
############### 2018314030 김단아 ################

################
#### R 코드 ####
################


rm(list=ls())
setwd("C:\\Users\\User\\Desktop\\18-2\\DM\\data") # harris.dat이 존재하는 위치로 지정해주어야한다.
getwd()


# make multiple linear regression function

regression.analysis <- function() {
  
  # import data file
  data = readline("Enter the data file name (with extension name) : " )
  cat("Select the data coding format(1 = 'a b c' or 2 = 'a,b,c'): ")
  fm = scan(n=1, quiet=TRUE)
  if(fm==1) {form = ""} else {form = ","}
  data = read.table(data, sep=form)
  
  # enter the Column number
  cat("Enter which column the response variable is recorded: ")
  num = scan(n=1, quiet=TRUE)
  
  # design matrix
  n = dim(data)[1]
  p = dim(data)[2]-1
  one = matrix(1, nrow=n, ncol=1)
  I = diag(n)
  
  y = as.matrix(data[,num])
  x = as.matrix(cbind(one,data[,-num]))
  H = x%*%solve(t(x)%*%x)%*%t(x)
  H0 = one%*%solve(t(one)%*%one)%*%t(one)
  
  # multiple regression result
  b = round(solve(t(x)%*%x)%*%t(x)%*%y, 3)
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
  write.table(b, outputname, sep= "   ", row.names=TRUE, col.names=FALSE, append=TRUE, quote=FALSE)
  
  cat("\n","ID, Actual values, Fitted values","\n","--------------------------------","\n",file = outputname,sep="", append=TRUE)
  write.table(y.values, outputname, sep= ", ", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
  
  cat("\n","Model Summary","\n","-------------","\n",file = outputname,sep="", append=TRUE)
  write.table(model.smy, outputname, sep= "", row.names=TRUE, col.names=FALSE, append=TRUE, quote=FALSE)
  
  cat("Output file has been saved in ",getwd(),"/",outputname,sep="")
}


regression.analysis()
harris.dat
1
1
HW2KimDA_R_output





#####################
#### Python 코드 ####
#####################



def regression_analysis():
  
  # import packages
  import numpy as np
  import numpy.linalg as lin
  import os as os
  import pandas as pd
  
  # prompt user to enter the data information
  data = input("Enter the data file name: ")
  fm = input("Select the data coding format(1 = 'a b c' or 2 = 'a,b,c'): ")
  num = int(input("Enter which column the response variable is recorded: "))-1
  
  if fm=='1': 
    form = ' '  
  else: 
    form= ','
  
  # import data
  data = pd.read_csv(data, header=None, sep=form)
  response = data[num]
  explanatory = data.drop(num, axis=1)
  
  # design matrix
  n = data.shape[0]
  p = data.shape[1]-1
  one = pd.DataFrame(np.ones((n,1)))
  I = np.eye(n)
  X = pd.concat([pd.DataFrame(one),explanatory], axis=1)
  Y = response
  H = X.dot(lin.inv(X.T.dot(X))).dot(X.T)
  H0 = one.dot(lin.inv(one.T.dot(one))).dot(one.T)
  inv = lin.inv(X.T.dot(X))
  b = inv.dot(X.T).dot(Y).round(3)
  
  # multiple regression result
  yhat = H.dot(Y).round(1)
  SST = Y.T.dot(I-H0).dot(Y)
  SSE = Y.T.dot(I-H).dot(Y)
  Rsquare = round(1 - SSE/SST, 4)
  MSE = round(SSE/(n-p), 4)
  
  # output file name
  outputname = input("Write the output file name you want to save (without extension name) : ")
  outputname = outputname+'.txt'
  
  # outport the result
  with open(outputname,"w") as text_file:
    
    print("Coefficients", file=text_file)
  print("-------------", file=text_file)
  for i in range(p+1):
    if i==0:
    print("Constant:", b[i],sep="   ", file=text_file)
  else: 
    print("Beta",i,":   ",b[i],sep="", file=text_file)
  print("",file=text_file)
  
  print("ID, Actual values, Fitted values", file=text_file)
  print("--------------------------------", file=text_file)
  for i in range(n):
    print(i+1, Y[i], yhat[i], sep=", ", file=text_file)
  print("",file=text_file)
  
  print("Model Summary", file=text_file)
  print("-------------", file=text_file)
  print("R-square = ", Rsquare, sep="", file=text_file)
  print("MSE = ", MSE, sep="", file=text_file)
  

