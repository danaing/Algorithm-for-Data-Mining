rm(list=ls())

# check working library
mylib = function() {
  mylib = readline('Write the location of the data file. : ')
  setwd(mylib)
  cat('Working directory is now', getwd(),'\n')
}

read = function(){
  data = readline("Enter the data file name (with extension name) : " )
  cat("Select the data coding format(1 = 'a b c' or 2 = 'a,b,c'): ")
  fm = scan(n=1, quiet=TRUE)
  if(fm == 1){form=""} else {form=","}
  read.table(data, sep=form, header = TRUE)
}

mylib()
C:\\Users\\User\\Desktop\\18-2\\DM\\data


# import training & test data file
cat('Import the dataset of TRAINING','\n')
train = read()
titanic.csv
2
train

test = read()
titanic.csv
2

# enter the Column number
cat("Enter which column the response variable is recorded: ")
num = scan(n=1, quiet=TRUE)

# nclass of response variable
k = length(unique(train[,num]))      # Assume that values of the class variable are integers starting with 1

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
GoS = matrix(ncol=2, nrow=p)
s = c()

for (j in 1:p) {
  # Xi variable
  value = sort(unique(x[,j]))
  c = length(value)
  combi = expand.grid(rep(list(0:1), c))
  t = nrow(combi)/2
  goodness <- c() 
  cat('Start Greedy Search for', colnames(train)[-num][j], 'variable.','\n')
  cat('Repeat ',t-1 ,'times.','\n')
  
  for (i in 1:(t-1)) {
    # 2^(c-1)-1 time repetition
    s = value[as.numeric(combi[i+1,])*c(1:c)]
    cat('Split at', as.character(s), '\t')
    # split
    t1 = subset(y, x[,j] %in% s) ; t2 = subset(y, !(x[,j] %in% s))
    n1 = length(t1) ; n2 = length(t2)
    
    #Gini impurity and Goodness of split
    imp_t1 <- 1-(sum(t1==classes[1])/n1)^2-(sum(t1==classes[2])/n1)^2
    imp_t2 <- 1-(sum(t2==classes[1])/n2)^2-(sum(t2==classes[2])/n2)^2
    imp_t <- 1-(sum(y==classes[1])/n)^2-(sum(y==classes[2])/n)^2
    goodness[i] <- imp_t - n1/n*imp_t1 - n2/n*imp_t2
    cat('Goodness of split is ', round(goodness[i], 4), '\n')
  }
  GoS[j,] = c(which.max(goodness), round(max(goodness),4))
  cat('\n')
}

# result of split point and spliting
result <- GoS
rownames(result) <- c(colnames(train)[-num])
colnames(result) <- c('Index','Goodness of Split')
print(result)
cat('\n')

j = which.max(result[,2])
var_1 = rownames(result)[j]
t = result[j,1]

value = sort(unique(x[,j]))
c = length(value)
combi = expand.grid(rep(list(0:1), c))
split_1 = as.character(value[as.numeric(combi[t+1,])*c(1:c)])

# training set
cat('Node 1 is ', var_1, 'and split at',split_1, '\t')
t1 = subset(y, x[,j] %in% split_1); t2 = subset(y, !(x[,j] %in% split_1))
n1 = length(t1) ; n2 = length(t2)
table(t1)
c1 = as.character(classes[which.max(table(t1))]) ; c2 = as.character(classes[which.max(table(t2))])
class = rep(0,n)
class[x[,j] %in% split_1] <- c1
class[!(x[,j] %in% split_1)] <- c2

# test set
class.t = rep(0,n.t)
class.t[x.t[,j] %in% split_1] <- c1
class.t[!(x.t[,j] %in% split_1)] <- c2

# Output setting
predict = cbind(c(1:n), as.character(y), class)
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
cat("\t", 'Node 1: ', var_1, ' in ', split_1, ' (', table(y)[1], ',', table(y)[2], ')', '\n', file = outputname, sep="", append=TRUE)
cat("\t", 'Node 2: ', c1, ' (', table(t1)[1], ',', table(t1)[2], ')', '\n', file = outputname, sep="", append=TRUE)
cat("\t", 'Node 3: ', c2, ' (', table(t2)[1], ',', table(t2)[2], ')', '\n','\n', file = outputname, sep="", append=TRUE)

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



numer = as.numeric(strsplit(readline('Numerical attributes :'),split=",")[[1]])