
rm(list=ls())
mylib = "C:\\Users\\User\\Desktop\\18-2\\DM\\data"  # veh.dat이 존재하는 위치로 지정해주어야한다.
setwd(mylib)
getwd()

# For 3D graphic output 
#install.packages("rgl")
library(rgl)



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
  row.names(model.smy) = c("R-square = ","MSE = ")
  
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
  
  # choose (i)LDA or (ii)QDA or (iii)RDA
  cat("Enter 1 to implement LDA, 2 to implement QDA or 3 to implement RDA")
  choice = scan(n=1, quiet=TRUE)
  
  # nclass of response variable
  k = length(unique(train[,num]))      # Assume that values of the class variable are integers starting with 1

  # set prior
  cat('It has',k,'Classes.','\n',
      'Input each priors with ascending order of Class index.','\n',
      'ex) Priors of 3 class size = 1/3, 1/3, 1/3','\n')
  prior = strsplit(readline('Priors :'),split=",")[[1]]
  
  # Error type 1
  if (length(prior)!= k) {cat('Your Prior input does not correspond with class size.','\n') 
    ; prior = paste(rep(1, k),'/',k,sep="")
    ; cat('So equal prior is given such as', prior)}
  
  prior = as.vector(sapply(prior, function(x){eval(parse(text=prior))}))
  # Error type 2
  if (length(sum(prior))!= 1) {cat('Sum of your prior input is not equal to 1','\n')
    ; prior = paste(rep(1, k),'/',k,sep="")
    ; cat('So equal prior is given such as', prior)}
  
  # Basic vectors
  n = nrow(train)
  n.t = nrow(test)
  x = t(train[-num])
  x.t = t((test[,-num]))
  y = train[,num]
  y.t = test[,num]
  p = ncol(train)-1
  
  means = t(as.matrix(aggregate(train[-num], train[num], mean)[,-1]))
  sk = lapply(lapply(split(train,train[,num]), function(x){x[,-num]}), cov)
  sp_nonsum = lapply( sk, function(x){ (nrow(x)-1)*x/(n-k) })
  sp = Reduce('+', sp_nonsum)

  # (i) LDA
  if(choice==1) {
    
    # Training data
    d = t(means)%*%solve(sp)%*%x + matrix(rep(diag(-(1/2)*t(means)%*%solve(sp)%*%means),each=n),ncol=n,byrow=TRUE) + matrix(rep(log(prior),each=n),ncol=n, byrow=TRUE)
    c = as.matrix(apply(d, 2, function(x) which(x==max(x))))
    # Test data
    d.t = t(means)%*%solve(sp)%*%x.t + matrix(rep(diag(-(1/2)*t(means)%*%solve(sp)%*%means),each=n.t),ncol=n.t,byrow=TRUE) + matrix(rep(log(prior),each=n.t),ncol=n.t, byrow=TRUE)
    c.t = as.matrix(apply(d, 2, function(x) which(x==max(x))))
    
    # (ii) QDA
  } else 
    if(choice==2) {

    # Training data
    d = matrix(nrow=k, ncol=n)
    for (i in 1:k) {
      d[i,] = rep(-(1/2)*log(det(sk[[i]])),n) + diag(-(1/2)*t(apply(x,2,function(x){x-means[,i]}))%*%solve(sk[[i]])%*%apply(x,2,function(x){x-means[,i]})) + rep(log(prior[i]),each=n)
    }
    c = as.matrix(apply(d, 2, function(x) which(x==max(x))))
    # Test data
    d.t = matrix(nrow=k, ncol=n.t)
    for (i in 1:k) {
      d.t[i,] = rep(-(1/2)*log(det(sk[[i]])),n.t) + diag(-(1/2)*t(apply(x.t,2,function(x){x-means[,i]}))%*%solve(sk[[i]])%*%apply(x.t,2,function(x){x-means[,i]})) + rep(log(prior[i]),each=n.t)
    }
    c.t = as.matrix(apply(d.t, 2, function(x) which(x==max(x))))
   
    # (iii) RDA
   } else 
     if(choice==3) {
    
    #Alpha & gamma selection by 0.05
    alpha = seq(0,1,0.05)
    gamma = seq(0,1,0.05)
    d_temp = matrix(nrow=k, ncol=n)
    s_rda = matrix(nrow=p, ncol=p)
    sigma = mean(diag(sp))
    result = c()

    for (a in alpha) {
      for (g in gamma) {
        for (i in 1:k) {
          s_rda = a*sk[[i]]+(1-a)*(g*sp+(1-g)*sigma*diag(p))
          d_temp[i,] = rep(-(1/2)*log(det(s_rda)),n) + diag(-(1/2)*t(apply(x,2,function(x){x-means[,i]}))%*%solve(s_rda)%*%apply(x,2,function(x){x-means[,i]})) + rep(log(prior[i]),each=n)
        }
        c_temp = as.matrix(apply(d_temp, 2, function(x) which(x==max(x))))
        #table_temp = table(y, c_temp)
        #accuracy_temp = sum(table_temp[row(table_temp)==col(table_temp)])/sum(table_temp)
        accuracy_temp = sum(c_temp==y)/n
        accuracy_rate = c(a, g, accuracy_temp)
        result <- rbind(result, accuracy_rate)
      }
    }
    
    # Optimal parameters
    optimal_loc = as.numeric(which(result[,3]==max(result[,3], na.rm = TRUE)))
    optimal = result[optimal_loc,]
    colnames(optimal) = c('alpha', 'gamma', 'accuracy rate')
    rownames(optimal) = paste(rep(1:length(optimal_loc)),'parameter')
    
    cat('Optimal parameters are given as follows.', '\n') ; print(optimal)
    
    # Training data
    a = optimal[1,1]
    g = optimal[1,2]
    d = matrix(nrow=k,ncol=n)
    s_rda_op = matrix(nrow=p, ncol=p)
    for (i in 1:k){
      s_rda_op = a*sk[[i]]+(1-a)*(g*sp+(1-g)*sigma*diag(p))
      d[i,] = rep(-(1/2)*log(det(s_rda_op)),n) + diag(-(1/2)*t(apply(x,2,function(x){x-means[,i]}))%*%solve(s_rda_op)%*%apply(x,2,function(x){x-means[,i]})) + rep(log(prior[i]),each=n)
    }
    c <- as.matrix(apply(d, 2, function(x) which(x==max(x))))
    
    # Test data
    d.t <- matrix(nrow=k, ncol=n.t)
    for (i in 1:k) {
      s_rda_op = a*sk[[i]]+(1-a)*(g*sp+(1-g)*sigma*diag(p))
      d.t[i,] = rep(-(1/2)*log(det(s_rda_op)),n.t) + diag(-(1/2)*t(apply(x.t, 2, function(x){x-means[,i]}))%*%solve(s_rda_op)%*%apply(x.t, 2, function(x){x-means[,i]})) + rep(log(prior[i]),each=n.t)
    }
    c.t <- as.matrix(apply(d.t, 2, function(x) which(x==max(x)))) 
    } else warning ('Choose 1 for LDA, 2 for QDA or 3 for RDA') 
  
  
  # Output setting

  cat('Please enter the number of the rows you can have.')
  out_num = as.numeric(readline('Enter the number : '))
  
  if (choice <= 3) {
    predict = head(cbind(c(1:n), y, c), n=out_num)
    table = table(y, c, dnn=c("Actual Class","Predicted Class"))
    accuracy = sum(diag(table))/sum(table)
    Predict.t = head(cbind(c(1:n.t),y.t,c.t),n=out_num)
    table.t = table(y.t, c.t, dnn=c("Actual Class","Predicted Class"))
    accuracy.t = sum(diag(table.t))/sum(table.t)
  
    # make output file
    outputname = readline("Write the output file name you want to save (without extension name) : ")
    outputname = paste(outputname,".txt",sep="")
    
    cat("ID, Actual class, Resub pred", "\n", "-----------------------------", "\n", file = outputname, sep="")
    write.table(predict, outputname, sep= ", ", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
    cat('(continue)','\n','\n', file = outputname, sep="", append = TRUE)
    cat('Confusion Matrix (Resubstitution)', "\n", "----------------------------------", "\n",file = outputname,sep="", append=TRUE)
    capture.output(print(table), file=outputname, append=TRUE)
    cat("\n", "Model Summary (Resubstitution)", "\n", "------------------------------", "\n",file = outputname, sep="", append=TRUE)
    cat("Overall accuracy: ", round(accuracy, 3), "\n\n",file = outputname, sep="", append=TRUE)
    cat("ID, Actual class, Test pred", "\n", "-----------------------------", "\n",file = outputname,sep="", append=TRUE)
    write.table(Predict.t, file=outputname, sep= ", ", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
    cat('(continue)',"\n",'\n', file = outputname, sep="", append = TRUE)
    cat('Confusion Matrix (Test)', "\n", "----------------------------------", "\n",file = outputname,sep="", append=TRUE)
    capture.output(print(table.t), file=outputname,append=TRUE)
    cat("\n", "Model Summary (Test)", "\n", "------------------------------", "\n",file = outputname,sep="", append=TRUE)
    cat("Overall accuracy: ", round(accuracy.t, 3), "\n" ,file = outputname,sep="", append=TRUE)
    cat('\n', "Chosen tuning parameters are\n alpha = ",a,"\n gamma = ",g,file=outputname, append=TRUE)
   
    # Plotting Test Accuracy
    if (choice==3) {
      library(rgl)
      X1 = rep(alpha,each=length(alpha))
      X2 = rep(gamma,length(gamma))
      Y  = as.vector(result[,3])
      plot3d(X1,X2,Y,type='s',size=1.5,col = 'red')
      bgplot3d({
        plot.new()
        title(main ="TEST ACCURACY" ,line = 2)
        mtext(side = 1,"X1 : Alpha,    X2 : Gamma,    Y : Test Accuracy", line = 3)
      })}
    
    cat("Output file has been saved in ",getwd(),"/",outputname,sep="")}  
  
  }



HW4 = function() {
  cat('Enter 1 to run Regression.','\n',
      'Enter 2 to run Classification.')
  ans = scan(n=1, quiet=TRUE)
  if(ans == 1){regression()} else {classification()}
}



HW4()
2
veh.dat
2
vehtest.dat
2
19
3
1/4,1/4,1/4,1/4
5
HW4KimDA_class_R_output


