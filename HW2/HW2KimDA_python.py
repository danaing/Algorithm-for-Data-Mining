
# coding: utf-8

# set the working directory
import os as os
os.chdir("C:/Users/User/Desktop/18-2/DM/data")
print(os.getcwd())


def regression_analysis():
    
    # import packages
    import numpy as np
    import numpy.linalg as lin
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