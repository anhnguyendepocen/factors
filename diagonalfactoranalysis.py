#Project    :  Resources, Diagonal Factor Analysis
#Description:  this .py file produces diagonal factor analysis for an arbitrary full rank
#              measurement system (the code crashes if the measurement system is not full rank)
#
#Basics:       three functions,                 (1) standardize data and obtain the correlation matrix 
#                                               (2) scree plot to determine the number of factors
#                                               (3) diagonal factor analysis
#                                               (*) option to avoid displaying scree plot
#                                               (*) option to manually input the number of factors
#              input,                           a linearly independent measurement system as an array
#              output,                          an array called fscores containning factor scores	
#This version: 06/19/2014
#This .py file:   Jorge Luis Garcia
#This project :   Jorge Luis Garcia

import numpy as np
import matplotlib.pyplot as mplot

#Correlation matrix (standardize data) 
def stdcorrdata(tofactordata):
    N, K = list(tofactordata.shape)
    N = float(N)
    K = int(K)
    #Standardize each measure
    for k in range(0, K):
        tofactordata[:,k] = (tofactordata[:,k] - np.mean(tofactordata[:,k]))/np.std(tofactordata[:,k])
    #Compute correlation matrix
    CM = np.dot(np.transpose(tofactordata),tofactordata)/(N)
    return CM
    
#Number of factors through Scree Plot rule (>1 eigen values of std correlation matrix)
def screeplot(tofactordata,display=False):
    N, K = list(tofactordata.shape)
    CM = stdcorrdata(tofactordata)
    #Compute eigenvalues, produce scree plot, define number of factors
    eigCM = np.linalg.eigvals(CM)
    
    #Plot
    if display==True:    
        mplot.plot([k for k in range(1,K+1)], np.sort(eigCM)[::-1], 'r-')
        mplot.plot([k for k in range(1,K+1)], np.ones(len(eigCM)), 'k-')
        #mplot.title("")
        #mplot.xlabel('Number of Factors');
        #mplot.ylabel('Eigenvalues');
        #mplot.savefig("ScreePlot_one_ccog.png")
        #mplot.close()    
    
    eigCM = np.sort(eigCM)[::-1]
    #Define and display number of factors
    factorn = max(sum(abs(eigCM) > 1),1)
    return factorn        
        
#Diagonal factor analysis
def diagonalfac(tofactordata,factorn=None):
    print str(tofactordata.shape[1]) + " number of measures"
    print str(np.linalg.matrix_rank(tofactordata)) + " are linearly independent"
    print str(tofactordata.shape[1]) + " should be identical to " + str(np.linalg.matrix_rank(tofactordata)) + " for the code to finish with success"
    if factorn is None:
        factorn = screeplot(tofactordata,display=False)
    N, K = list(tofactordata.shape)
    #Compute standardized correlation matrix
    CM = stdcorrdata(tofactordata)
    #Create a matrix storing the factor loadings and factor loadings matrices
        
    floads  = np.zeros([K,factorn])
    floadsm = np.zeros([K,K,factorn])
    
    #Ingredientes to extraxt factors
    for j in range (0,factorn):
        #sum of squares of the correlation matrix    
        ssCM = (CM**2).sum(axis=0)
        #Index variable with maximum correlation
        maxcorr2 = np.argmax(ssCM)
        #Calculate and store factor loadings
        floads[:,j] = CM[:,maxcorr2]
        #Calculate and store residual matrix
        floadsm[:,:,j] = np.dot((CM[:,maxcorr2]).reshape(-1,1),(CM[:,maxcorr2]).reshape(1,-1))
        #Residualize the correlation matrix 
        CM -= floadsm[:,:,j]
    print floads
    #Define and return factors
    fscores = np.dot(np.dot(tofactordata,floads),np.linalg.inv(np.dot(np.transpose(floads),floads)))
    return fscores