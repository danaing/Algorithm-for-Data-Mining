
rm(list=ls())
mylib = "C:\\Users\\User\\Desktop\\18-2\\DM\\data"  # veh.dat이 존재하는 위치로 지정해주어야한다.
setwd(mylib)
getwd()



### importing dataset function ###

read = function(){
  data = readline("Enter the data file name (with extension name) : " )
  cat("Select the data coding format(1 = 'a b c' or 2 = 'a,b,c'): ")
  fm = scan(n=1, quiet=TRUE)
  if(fm == 1){form=""} else {form=","}
  read.table(data, sep=form)
}



### regressoin function ###

regression <- function() {
  
  # import data file
  data = read()
  
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



### classification function ###

classification <- function() {
  
  # import training & test data file
  cat('Import the dataset of TRAINING','\n')
  train = read()
  cat('Import the dataset of TEST','\n')
  test = read()
  
  # enter the Column number
  cat("Enter which column the response variable is recorded: ")
  num = scan(n=1, quiet=TRUE)
  
  # choose (i)LDA or (ii)QDA
  cat("Enter 1 to implement LDA or enter 2 to implement QDA")
  choice = scan(n=1, quiet=TRUE)
  
  # set prior
  nclass = length(unique(train[,num]))      # Assume that values of the class variable are integers starting with 1
  
  cat('It has',nclass,'Classes.','\n',
      'Input each priors with ascending order of Class index.','\n',
      'ex) Priors of 3 class size = 1/3, 1/3, 1/3','\n')
  prior = strsplit(readline('Priors :'),split=",")[[1]]
  
  if (length(prior)!= nclass) {cat('Your Prior input does not correspond with class size.','\n') 
    ; prior = paste(rep(1, nclass),'/',nclass,sep="")
    ; cat('So equal prior is given such as', prior)}
  
  # Basic vectors
  N = nrow(train)
  n = nrow(test)
  x = t(train[-num])
  x.test = t(as.matrix(test[,-num]))
  prior = as.matrix(sapply(prior,function(x){eval(parse(text=x))}))
  means = t(as.matrix(aggregate(train[-num], train[num], mean)[,-1]))
  cov = lapply(lapply(split(train,train[,num]), function(x){x[,-num]}), cov)
  covs = lapply( cov, function(x){ (nrow(x)-1)*x /(N - nclass) })
  sp = Reduce('+', covs)
  
  # (i) LDA
  if(choice==1) {
    d.resub = matrix(0,nrow=N,ncol=nclass)
    d.test = matrix(0,nrow=n,ncol=nclass)
    
    for(i in 1:nclass) {
      t0 = t(t(means[,i])%*%solve(sp)%*%x)
      t1 = t(t(means[,i])%*%solve(sp)%*%x.test)
      t2 = t((-0.5*t(means[,i])%*%solve(sp)%*%means[,i] + log(prior[i]))%*%matrix(1,nrow=1,ncol=N))
      t3 = t((-0.5*t(means[,i])%*%solve(sp)%*%means[,i] + log(prior[i]))%*%matrix(1,nrow=1,ncol=n))
      d.resub[,i]=t0+t2
      d.test[,i]=t1+t3 }
    
    result.resub = cbind(train[,num], max.col(d.resub))
    result.test = cbind(test[,num], max.col(d.test))
    
  } else if(choice==2){
    
    # (ii) QDA
    d.resub = matrix(0,nrow=N,ncol=nclass)
    d.test = matrix(0,nrow=n,ncol=nclass)
    
    for(i in 1:nclass) {
      t0 = -0.5*log(det(cov[[i]]))
      for(j in 1:N) {
        t1 = -0.5*t(x[,j]-means[,i])%*%solve(cov[[i]])%*%(x[,j]-means[,i]) + log(prior[i])
        d.resub[j,i]=t0+t1 }
      for(j in 1:n) {
        t2 = -0.5*t(x.test[,j]-means[,i])%*%solve(cov[[i]])%*%(x.test[,j]-means[,i]) + log(prior[i])
        d.test[j,i]=t0+t2 }
    }
    
    result.resub = cbind(train[,num], max.col(d.resub))
    result.test = cbind(test[,num], max.col(d.test))  
  } else warning ('Choose 1 for LDA or 2 for QDA')
  
  # Confusion Matrix
  t1 = table(result.resub[,1], result.resub[,2], dnn=c("Actual Class","Predicted Class"))
  t2 = table(result.test[,1], result.test[,2], dnn=c("Actual Class","Predicted Class"))
  
  # Accuracy
  accuracy.resub = sum(result.resub[,1] == result.resub[,2])/N
  accuracy.test = sum(result.test[,1] == result.test[,2])/n
  
  # make output file
  outputname = readline("Write the output file name you want to save (without extension name) : ")
  outputname = paste(outputname,".txt",sep="")
  
  
  cat("ID, Actual class, Resub pred", "\n", "--------------------------", "\n", file = outputname, sep="")
  write.table(cbind(c(1:N), result.resub), outputname, sep= ", ", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
  
  cat("\n", "Confusion Matrix (Resubstitution)", "\n", "--------------------------------", "\n",file = outputname,sep="", append=TRUE)
  capture.output(print(t1),file=outputname,append=TRUE)
  
  cat("\n", "Model Summary (Resubstitution)", "\n", "--------------------------------", "\n",file = outputname,sep="", append=TRUE)
  cat("Overall accuracy: ", round(accuracy.resub,3), "\n\n",file = outputname, sep="", append=TRUE)
  
  
  cat("ID, Actual class, Test pred", "\n", "--------------------------", "\n",file = outputname,sep="", append=TRUE)
  write.table(cbind(c(1:n), result.test), file = outputname, sep= ", ", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
  
  cat("\n", "Confusion Matrix (Test)", "\n", "--------------------------------", "\n",file = outputname,sep="", append=TRUE)
  capture.output(print(t2),file=outputname,append=TRUE)
  
  cat("\n", "Model Summary (Test)", "\n", "--------------------------------", "\n",file = outputname,sep="", append=TRUE)
  cat("Overall accuracy: ", round(accuracy.test,3), "\n" ,file = outputname,sep="", append=TRUE)
  
  cat("Output file has been saved in ",getwd(),"/",outputname,sep="")
}



HW3 = function() {
  cat('Enter 1 to run Regression.','\n',
      'Enter 2 to run Classification.')
  ans = scan(n=1, quiet=TRUE)
  if(ans == 1){regression()} else {classification()}
}





# for 1.b
HW3()
1
harris.dat
1
1
HW3KimDA_Reg_R_output


# for classification of veh.dat and vehtest.dat  
HW3()
2    # classification
veh.dat    # training dataset 
2    # form
vehtest.dat    # test dataset
2    # form
19    # response var
2    # QDA
1/4,1/4,1/4,1/4    # arbitrary set equal prior
HW3KimDA_R_output    # ouput file name
