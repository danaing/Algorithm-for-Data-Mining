# Generate a training sample of size 1000
# train set 1000
set.seed(1234)
a1 = matrix(c(2, 2),nrow=2)
a2 = matrix(c(3, -3),nrow=2)
X1 = rnorm(1000, 0, 1)
X2 = rnorm(1000, 0, 1)
X = rbind(X1, X2)
Z = rnorm(1000, 0, 1)
sigmoid = function(x){1/(1+exp(-x))}
Y = c(sigmoid(t(a1)%*%X) + (t(a2)%*%X)^2 + 0.3*Z)
train = data.frame(X1,X2,Y)
head(train)

# test set 1000
set.seed(5678)
a1 = matrix(c(2, 2),nrow=2)
a2 = matrix(c(3, -3),nrow=2)
X1 = rnorm(1000, 0, 1)
X2 = rnorm(1000, 0, 1)
X = rbind(X1, X2)
Z = rnorm(1000, 0, 1)
sigmoid = function(x){1/(1+exp(-x))}
Y = c(sigmoid(t(a1)%*%X) + (t(a2)%*%X)^2 + 0.3*Z)
test = data.frame(X1,X2,Y)
head(test)

## a. Plot the surface of responses
#install.packages('rgl')
require(rgl)

X1 = seq(from=-3, to=3, length.out = 100)
X2 = seq(from=-3, to=3, length.out = 100)
a1 = matrix(c(2, 2),nrow=2)
a2 = matrix(c(3, -3),nrow=2)

Y = function(x,y) {
  X = rbind(x, y)
  sigmoid(t(a1)%*%X) + (t(a2)%*%X)^2 
}

Z <- outer(X1, X2, Y)

persp3d(X1, X2, Z, theta=50, phi=25, expand=0.75, col='blue',
        ticktype="detailed", xlab="X1", ylab="X2", zlab="Response",axes=TRUE)
surface3d(X1, X2, Z, back = "lines")
surface3d(X1, X2, Z, front = "lines")


## Neural network analysis using nnet 
#install.packages('nnet')
library(nnet)
#install.packages('devtools')
library(devtools)
source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')
library(reshape2)

# hidden nodes from 1 to 10
train_error = c()
test_error = c()

for (i in 1:10) {
  nn <- nnet(Y~., data=train, size=i, linout = TRUE)
  print(summary(nn))
  
  #### validation
  # train
  pred <- predict(nn, train[,-3])
  train_error[i] = sum((train$Y-pred)^2)/length(pred)
  
  # test
  pred_t <- predict(nn, test[,-3])
  test_error[i] = sum((test$Y-pred_t)^2)/length(pred_t)
}

## plot the error curves 
plot(train_error, type = 'o', col='blue', ylab='Error', xlab='the number of hidden node')
lines(test_error, col='red', type='o', lty = 2)
legend('topright', c('train set', 'test set'), col = c('blue', 'red'), lty = c(1,2))
