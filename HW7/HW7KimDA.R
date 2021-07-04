###############################################
##### 0. Checking the working environment #####
###############################################

rm(list=ls())

# check working library
mylib = function() {
  mylib = readline('Write the location of the data file. : ')
  setwd(mylib)
  cat('Working directory is now', getwd(),'\n')
}

# check installed package
is.install = function(package) {
  if(!is.element(package, installed.packages()[,1])) {install.packages(package)}
  else {cat(package,"is already installed. \n")} 
}


######################################
##### importing dataset function #####
######################################

read = function(){
  data = readline("Enter the data file name (with extension name) : " )
  cat("Select the data coding format(1 = 'a b c' or 2 = 'a,b,c'): ")
  fm = scan(n=1, quiet=TRUE)
  if(fm == 1){form=""} else {form=","}
  read.table(data, sep=form)
}

##################################
##### 1. regressoin function #####
##################################

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




######################################
##### 2. classification function #####
######################################

classification <- function() {
  
  ###################### General background #######################
  
  # import training & test data file
  cat('Import the dataset of TRAINING','\n')
  train = read()
  cat('Import the dataset of TEST','\n')
  test = read()
  
  # enter the Column number
  cat("Enter which column the response variable is recorded: ")
  num = scan(n=1, quiet=TRUE)
  
  # nclass of response variable
  k = length(unique(train[,num]))      # Assume that values of the class variable are integers starting with 1
  
  # choose (i)LDA or (ii)QDA or (iii)RDA (iv)Logistic Regression (v)Naive Bayes (vi)1-level decision tree
  cat("Enter 1 for LDA, 2 for QDA, 3 for RDA, 4 for Logistic Regression, 5 for Naive Bayes, 6 for 1-level decision tree.")
  choice = scan(n=1, quiet=TRUE)
  
  # do not prompt (iv) (v) (vi) when the data has more than 2 classes.
  if(k>2&choice==4|k>2&choice==5|k>2&choice==6) {stop('Cannot run chosen method since data has more than 2 classes.') }
  

  ################### Analysing data by chosen method ###################
  
  # Prepare for (i)LDA (ii)QDA (iii)RDA
  if ((choice==1)|(choice==2)|(choice==3)) {
    
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
  }
  
  ###############
  ### (i) LDA ###
  ###############
  if(choice==1) {
    # Training data
    d = t(means)%*%solve(sp)%*%x + matrix(rep(diag(-(1/2)*t(means)%*%solve(sp)%*%means),each=n),ncol=n,byrow=TRUE) + matrix(rep(log(prior),each=n),ncol=n, byrow=TRUE)
    c = as.matrix(apply(d, 2, function(x) which(x==max(x))))
    # Test data
    d.t = t(means)%*%solve(sp)%*%x.t + matrix(rep(diag(-(1/2)*t(means)%*%solve(sp)%*%means),each=n.t),ncol=n.t,byrow=TRUE) + matrix(rep(log(prior),each=n.t),ncol=n.t, byrow=TRUE)
    c.t = as.matrix(apply(d, 2, function(x) which(x==max(x))))
  }  
  
  ################
  ### (ii) QDA ###
  ################
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
  } 
  
  #################
  ### (iii) RDA ###
  ################# 
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
  } 
  
  #########################################
  ### Output setting for (i) (ii) (iii) ###
  #########################################
  
  if (choice <= 3) {
    
    cat('Please enter the row number you want in the output file.')
    out_num = as.numeric(readline('Enter the number : '))
    
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
    if (choice == 3) {
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
    
    cat("Output file has been saved in ",getwd(),"/",outputname,sep="") 
    }  
  
  ################################
  ### (iv) Logistic Regression ###
  ################################
  
  if (choice == 4) {
    # Basic matrix for Logistic
    n = nrow(train)
    n.t = nrow(test)
    x = cbind(rep(1,n), as.matrix(train[-num]))
    x.t = cbind(rep(1,n.t), as.matrix(test[-num]))
    y = train[,num]
    y.t = test[,num]
    p = ncol(train)
    
    library(maxLik)
    B = matrix(0, nrow=(p))
    loglik <- function(B) { t(y-1)%*%x%*%B - sum(log(1+exp(x%*%B))) }
    Lik <- maxLik(loglik, start=matrix(0, ncol=1, nrow=p))
    b <- coef(Lik)
    cat('The result that maximizes Log-Likelihood of Logistic Regression(MLE method).\n')
    print(summary(Lik))
    
    cat('Enter the cutoff value from 0 to 1.')
    cutoff = scan(n=1, quiet=TRUE)
    
    #Training data
    prob <- round(exp(x%*%b)/(1+exp(x%*%b)), 3)
    class <- c()
    for(i in 1:n){if (prob[i]>=cutoff) {class[i] = 2} else {class[i] = 1}}
    
    #Test data
    prob.t <- round(exp(x.t%*%b)/(1+exp(x.t%*%b)), 3)
    class.t <- c()
    for (i in 1:n.t) {if(prob.t[i]>=cutoff) {class.t[i] = 2} else {class.t[i] = 1}}
    
    # Output setting
    predict = cbind(c(1:n), y, class, prob)
    table = table(y, class, dnn=c("Actual Class","Predicted Class"))
    accuracy = sum(diag(table))/sum(table)
    sensi = table[2,2]/sum(table[2,])
    speci = table[1,1]/sum(table[1,])
    
    Predict.t = cbind(c(1:n.t), y.t, class.t, prob.t)
    table.t = table(y.t, class.t, dnn=c("Actual Class","Predicted Class"))
    accuracy.t = sum(diag(table.t))/sum(table.t)
    sensi.t = table.t[2,2]/sum(table.t[2,])
    speci.t = table.t[1,1]/sum(table.t[1,])
    
    # make output file
    out_num = as.numeric(readline('Please enter the maximum output row you want to have in the output file. :' ))
    
    outputname = readline("Write the output file name you want to save (without extension name) : ")
    outputname = paste(outputname,".txt",sep="")
    
    cat("ID, Actual class, Resub pred, Pred Prob", "\n", "-----------------------------", "\n", file = outputname, sep="")
    write.table(head(predict, out_num), outputname, sep= ", ", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
    cat('(continue)','\n','\n', file = outputname, sep="", append = TRUE)
    cat('Confusion Matrix (Resubstitution)', "\n", "----------------------------------", "\n",file = outputname,sep="", append=TRUE)
    capture.output(print(table), file=outputname, append=TRUE)
    cat("\n", "Model Summary (Resubstitution)", "\n", "------------------------------", "\n",file = outputname, sep="", append=TRUE)
    cat("Overall accuracy = ", round(accuracy, 3), "\n",file = outputname, sep="", append=TRUE)
    cat("Sensitivity = ", round(sensi, 3), "\n",file = outputname, sep="", append=TRUE)
    cat("Specificity =", round(speci, 3), "\n\n",file = outputname, sep="", append=TRUE)
    
    cat("ID, Actual class, Test pred, Pred Prob", "\n", "-----------------------------", "\n",file = outputname,sep="", append=TRUE)
    write.table(head(Predict.t, out_num), file=outputname, sep= ", ", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
    cat('(continue)',"\n",'\n', file = outputname, sep="", append = TRUE)
    cat('Confusion Matrix (Test)', "\n", "----------------------------------", "\n",file = outputname,sep="", append=TRUE)
    capture.output(print(table.t), file=outputname,append=TRUE)
    cat("\n", "Model Summary (Test)", "\n", "------------------------------", "\n",file = outputname,sep="", append=TRUE)
    cat("Overall accuracy = ", round(accuracy.t, 3), "\n" ,file = outputname,sep="", append=TRUE)
    cat("Sensitivity = ", round(sensi.t, 3), "\n",file = outputname, sep="", append=TRUE)
    cat("Specificity =", round(speci.t, 3), "\n\n",file = outputname, sep="", append=TRUE)
    cat("Output file has been successfully saved in ",getwd(),"/",outputname,sep="") 
    } 
  
  #######################
  ### (v) Naive Bayes ###
  #######################
  
  if (choice == 5) {

    # set variable attributes
    cat('It has',k,'variables.','\n',
        'Input variable number which has numerical attribute.\n',
        'ex) 1, 3, 5 \n')
    numer = as.numeric(strsplit(readline('Numerical attributes :'),split=",")[[1]])
    
    # Basic vectors
    n = nrow(train)
    x = train[-num]
    y = train[,num]
    p = ncol(train)
    n.t = nrow(test)
    x.t = test[,-num]
    y.t = test[,num]

    # probability of training set
    t1 = train[train[,num]==1,]
    t2 = train[train[,num]==2,]
    prob_1 = matrix(1, nrow=n, ncol=p)
    prob_2 = matrix(1, nrow=n, ncol=p)
    
    for(i in 1:n) {
      for(j in c(1:p)[-num]) {
        if (train[i,j]=='?') {} else { 
            if (any(j==numer)) {
              prob_1[i,j] <- dnorm(train[i,j], mean(t1[,j]), sd(t1[,j]))
              prob_2[i,j] <- dnorm(train[i,j], mean(t2[,j]), sd(t2[,j]))
              } else { 
                t <- prop.table(table(y, x[[j]]), 1)
                prob_1[i,j] <- t[1,][colnames(t)==train[i,j]]
                prob_2[i,j] <- t[2,][colnames(t)==train[i,j]]
              }}}}
    
    prob_temp <- cbind(apply(prob_1,1,prod),apply(prob_2,1,prod))
    prob <- round(prob_temp/rowSums(prob_temp),3)
    class <- as.matrix(apply(prob, 1, function(x) which(x==max(x))))

    # probability of test set
    prob_1.t = matrix(1, nrow=n.t, ncol=p)
    prob_2.t = matrix(1, nrow=n.t, ncol=p)
    
    for(i in 1:n.t) {
      for(j in c(1:p)[-num]) {
        if (test[i,j]=='?') {} else { 
          if (any(j==numer)) {
            prob_1.t[i,j] <- dnorm(test[i,j], mean(t1[,j]), sd(t1[,j]))
            prob_2.t[i,j] <- dnorm(test[i,j], mean(t2[,j]), sd(t2[,j]))
          } else { 
            t <- prop.table(table(y, x[[j]]), 1)
            prob_1.t[i,j] <- t[1,][colnames(t)==test[i,j]]
            prob_2.t[i,j] <- t[2,][colnames(t)==test[i,j]]
          }}}}
    
    prob_temp.t <- cbind(apply(prob_1,1,prod),apply(prob_2,1,prod))
    prob.t <- round(prob_temp.t/rowSums(prob_temp.t),3)
    class.t <- as.matrix(apply(prob.t, 1, function(x) which(x==max(x))))
 
    # Output setting
    predict = cbind(c(1:n), y, class, apply(prob, 1, max))
    table = table(y, class, dnn=c("Actual Class","Predicted Class"))
    accuracy = sum(diag(table))/sum(table)
    sensi = table[2,2]/sum(table[2,])
    speci = table[1,1]/sum(table[1,])
    
    Predict.t = cbind(c(1:n.t), y.t, class.t, apply(prob.t, 1, max))
    table.t = table(y.t, class.t, dnn=c("Actual Class","Predicted Class"))
    accuracy.t = sum(diag(table.t))/sum(table.t)
    sensi.t = table.t[2,2]/sum(table.t[2,])
    speci.t = table.t[1,1]/sum(table.t[1,])

    # make output file
    out_num = as.numeric(readline('Please enter the maximum output row you want to have in the output file. :' ))
    
    outputname = readline("Write the output file name you want to save (without extension name) : ")
    outputname = paste(outputname,".txt",sep="")
    
    cat("ID, Actual class, Resub pred, Pred Prob", "\n", "-----------------------------", "\n", file = outputname, sep="")
    write.table(head(predict, out_num), outputname, sep= ", ", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
    cat('(continue)','\n','\n', file = outputname, sep="", append = TRUE)
    cat('Confusion Matrix (Resubstitution)', "\n", "----------------------------------", "\n",file = outputname,sep="", append=TRUE)
    capture.output(print(table), file=outputname, append=TRUE)
    cat("\n", "Model Summary (Resubstitution)", "\n", "------------------------------", "\n",file = outputname, sep="", append=TRUE)
    cat("Overall accuracy = ", round(accuracy, 3), "\n",file = outputname, sep="", append=TRUE)
    cat("Sensitivity = ", round(sensi, 3), "\n",file = outputname, sep="", append=TRUE)
    cat("Specificity =", round(speci, 3), "\n\n",file = outputname, sep="", append=TRUE)
    
    cat("ID, Actual class, Test pred, Pred Prob", "\n", "-----------------------------", "\n",file = outputname,sep="", append=TRUE)
    write.table(head(Predict.t, out_num), file=outputname, sep= ", ", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
    cat('(continue)',"\n",'\n', file = outputname, sep="", append = TRUE)
    cat('Confusion Matrix (Test)', "\n", "----------------------------------", "\n",file = outputname,sep="", append=TRUE)
    capture.output(print(table.t), file=outputname,append=TRUE)
    cat("\n", "Model Summary (Test)", "\n", "------------------------------", "\n",file = outputname,sep="", append=TRUE)
    cat("Overall accuracy = ", round(accuracy.t, 3), "\n" ,file = outputname,sep="", append=TRUE)
    cat("Sensitivity = ", round(sensi.t, 3), "\n",file = outputname, sep="", append=TRUE)
    cat("Specificity =", round(speci.t, 3), "\n\n",file = outputname, sep="", append=TRUE)
    cat("Output file has been successfully saved in ",getwd(),"/",outputname,sep="") 
    }
    

    ##################################
    ### (vi) 1-level decision tree ###
    ##################################
    
  if (choice == 6) {
      
    # Basic vectors
    classes = sort(unique(train[,num]))
    n = nrow(train)
    n.t = nrow(test)
    x = train[-num]
    x.t = test[,-num]
    y = train[,num]
    y.t = test[,num]
    p = ncol(train)-1
  
    # Greedy search for continuous variable
    GoS <- c()
    split <- c()
    s=0
    
    for (j in 1:p) {
      # Xi variable
      value=sort(unique(x[,j]))
      t=length(value)
      goodness <- c() 
      cat('Start Greedy Search for', colnames(train)[-num][j], 'variable.','\n')
      cat('Repeat ',t-1,'times.','\n')
      
      for (i in 1:(t-1)) {
        # t-1 time repetition
        s[i] = mean(value[c(i,(i+1))])
        
        cat('Split at',s[i], '\t')
        index_t1 = which(x[,j]<=s[i])
        index_t2 = which(x[,j]>s[i])
        
        t1 = y[index_t1] ; t2 = y[index_t2]
        n1 = length(y[index_t1]) ; n2 = length(y[index_t2])
        
        #Gini impurity
        imp_t1 <- 1-(sum(t1==classes[1])/n1)^2-(sum(t1==classes[2])/n1)^2
        imp_t2 <- 1-(sum(t2==classes[1])/n2)^2-(sum(t2==classes[2])/n2)^2
        
        #Goodness of split
        imp_t <- 1-(sum(y==classes[1])/n)^2-(sum(y==classes[2])/n)^2
        goodness[i] <- imp_t - n1/n*imp_t1 - n2/n*imp_t2
        cat('Goodness of split is ', round(goodness[i], 4), '\n')
      }
      GoS[j] = max(goodness)
      split[j] = s[which.max(goodness)]
      cat('\n')
    }
    
    # result of split point and spliting
    result = cbind(split, round(GoS,4))
    rownames(result) <- c(colnames(train)[-num])
    colnames(result) <- c('split','Goodness of Split')
    print(result)
    cat('\n')
    
    j=which.max(result[,2])
    var_1 = rownames(result)[j]
    split_1 = result[1,j]
    cat('Node 1 is ', var_1, 'and split at',split_1, '\t')
    index_t1 = which(x[,j]<=split_1)
    index_t2 = which(x[,j]>split_1)
    class = rep(0,n)
    class[index_t1]<-classes[1]
    class[index_t2]<-classes[2]
    t1 = y[index_t1] ; t2 = y[index_t2]
    
    # test set
    index_t1.t = which(x.t[,j]<=split_1)
    index_t2.t = which(x.t[,j]>split_1)
    class.t = rep(0,n.t)
    class.t[index_t1.t]<-classes[1]
    class.t[index_t2.t]<-classes[2]
    
    # Output setting
    predict = cbind(c(1:n), y, class)
    table = table(y, class, dnn=c("Actual Class","Predicted Class"))
    accuracy = sum(diag(table))/sum(table)
    sensi = table[2,2]/sum(table[2,])
    speci = table[1,1]/sum(table[1,])
    
    Predict.t = cbind(c(1:n.t), y.t, class.t)
    table.t = table(y.t, class.t, dnn=c("Actual Class","Predicted Class"))
    accuracy.t = sum(diag(table.t))/sum(table.t)
    sensi.t = table.t[2,2]/sum(table.t[2,])
    speci.t = table.t[1,1]/sum(table.t[1,])
    
    # make output file
    out_num = as.numeric(readline('Please enter the maximum output row you want to have in the output file. :' ))
    
    outputname = readline("Write the output file name you want to save (without extension name) : ")
    outputname = paste(outputname,".txt",sep="")
    
    cat("Tree Structure", "\n", file = outputname, sep="")
    cat("\t", 'Node 1: ', var_1, ' <= ', split_1, ' (', table(y)[1], ',', table(y)[2], ')', '\n', file = outputname, sep="", append=TRUE)
    cat("\t", 'Node 2: ', classes[1], ' (', table(t1)[1], ',', table(t1)[2], ')', '\n', file = outputname, sep="", append=TRUE)
    cat("\t", 'Node 3: ', classes[2], ' (', table(t2)[1], ',', table(t2)[2], ')', '\n','\n', file = outputname, sep="", append=TRUE)
    
    cat("ID, Actual class, Resub pred", "\n", "-----------------------------", "\n", file = outputname, sep="", append=TRUE)
    write.table(head(predict, out_num), outputname, sep= ", ", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
    cat('(continue)','\n','\n', file = outputname, sep="", append = TRUE)
    cat('Confusion Matrix (Resubstitution)', "\n", "----------------------------------", "\n",file = outputname,sep="", append=TRUE)
    capture.output(print(table), file=outputname, append=TRUE)
    cat("\n", "Model Summary (Resubstitution)", "\n", "------------------------------", "\n",file = outputname, sep="", append=TRUE)
    cat("Overall accuracy = ", round(accuracy, 3), "\n",file = outputname, sep="", append=TRUE)
    cat("Sensitivity = ", round(sensi, 3), "\n",file = outputname, sep="", append=TRUE)
    cat("Specificity = ", round(speci, 3), "\n\n",file = outputname, sep="", append=TRUE)
    
    cat("ID, Actual class, Test pred", "\n", "-----------------------------", "\n",file = outputname,sep="", append=TRUE)
    write.table(head(Predict.t, out_num), file=outputname, sep= ", ", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
    cat('(continue)',"\n",'\n', file = outputname, sep="", append = TRUE)
    cat('Confusion Matrix (Test)', "\n", "----------------------------------", "\n",file = outputname,sep="", append=TRUE)
    capture.output(print(table.t), file=outputname,append=TRUE)
    cat("\n", "Model Summary (Test)", "\n", "------------------------------", "\n",file = outputname,sep="", append=TRUE)
    cat("Overall accuracy = ", round(accuracy.t, 3), "\n" ,file = outputname,sep="", append=TRUE)
    cat("Sensitivity = ", round(sensi.t, 3), "\n",file = outputname, sep="", append=TRUE)
    cat("Specificity = ", round(speci.t, 3), "\n\n",file = outputname, sep="", append=TRUE)
    cat("Output file has been successfully saved in ",getwd(),"/",outputname,sep="") 
  }  
}



######################################
##### 3. Combining the functions #####
######################################

HW7 = function(){
  cat('Checking the working environment. \n')
  mylib()
  cat('Checking the packages required. \n')
  is.install('rgl')
  is.install('maxLik')
  ans = readline('Enter 1 to use Regression or 2 to use Classification : ')
  if (ans==1) { regression() }
  if (ans==2) {classification()}
}


####################################################
##### Implement Logistic Regression to pid.dat #####
####################################################

HW7()
C:\\Users\\User\\Desktop\\18-2\\DM\\data
2
harris.dat
1
harris.dat
1
5
6
5
HW7KimDA_tree_R_output


