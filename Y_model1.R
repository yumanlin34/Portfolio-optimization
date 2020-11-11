# Project 1 for CFRM 507
# Yuman Lin
# model 1 with Yeats

# model 1 (need to change for each model)
target1 = c(0.3, 0.09, 0.07, 0.2, 0.08, 0.06, 0.03, 0.03, 0.02, 0.02, 0.01, 0.02, 0.07,0,0)
target2 = c(0.24, 0.09, 0.06, 0.11, 0.05, 0.05, 0.03, 0.03, 0.01, 0.01, 0.01, 0.03, 0.26, 0, 0.02)
target3 = c(0.16, 0.06, 0.04, 0.08, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01, 0.02, 0.04, 0.35, 0.04, 0.08)

# constraint matrix
# weight vector
a = rep(1,15) #vector
a = t(a)

# first primary characteristic
# Pchar1 is the same as Pchar 3
Pchar1 = c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,0)
Pchar1 = t(Pchar1)
Pchar3 = c(0,0,0,0,0,0,0,1,1,1,0,0,0,0,1)
Pchar3 = t(Pchar3)

# for Amat, 3*55 with value 0
zero_vector=rep(0,55)
Amat1=rbind(a, Pchar1, Pchar3)                   # vector for non-zero
Amat2=rbind(zero_vector,zero_vector,zero_vector) # vector for zero
Amat_first=cbind(Amat1, Amat2)                        # vector for first part

# secondary characteristics
Schar1 = c(1,1,1,0,0,0.55,0.6,0,0,0,0,0,0,0,0)
Schar2 = c(0,0,1,0,0,0.1,1,0,0,0,0,0,0,0,0)
Schar3 = c(0,0,0,0,0,0,0,1,1,1,0,0,0,0,0)
Schar4 = c(0,0,0,0,0,0,0,0,0,0,1,1,0,0,0)
Schar5 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0)
Schar1 = t(Schar1)
Schar2 = t(Schar2)
Schar3 = t(Schar3)
Schar4 = t(Schar4)
Schar5 = t(Schar5)

# for d+ 5*5
dplus = diag(5)

# for d- 5*5
dminus = -diag(5)

Amat4 = rbind(Schar1, Schar2, Schar3, Schar4, Schar5) # vector for non-zero

# 5*45 with value 0
zero_vector1 = rep(0,45)
Amat5 = rbind(zero_vector1,zero_vector1,zero_vector1,zero_vector1,zero_vector1)

Amat_second = cbind(Amat4, dplus, dminus, Amat5) # vector for second part

# for home office model target, deltap and deltam, 15*15 identical matrix
deltap = diag(15)
deltam = -diag(15)

zero_matrix1 = matrix(0, nrow=15, ncol=10, byrow=TRUE) # 15*10 zero matrix
zero_matrix2 = matrix(0, nrow=15, ncol=15, byrow=TRUE) # 15*15 zero matrix

Amat_third = cbind(deltap, zero_matrix1, deltap, deltam, zero_matrix2)

# the fourth part for the relationship between weight and binary variable
zero_matrix3 = matrix(0, nrow=15, ncol=40, byrow=TRUE) # 15*10 zero matrix
Amat_fourth = cbind(deltap, zero_matrix3, deltam)

# last part: for the number of binary variable = 1
zero_matrix4 = matrix(0, nrow=1, ncol=55, byrow=TRUE)
Amat_fifth = cbind(zero_matrix4, a)

#====================================================Amat: 39*70==================================================
Amat = rbind(Amat_first, Amat_second, Amat_third, Amat_fourth, Amat_fifth)

# change the col and row name
rownames(Amat) <- c(1:dim(Amat)[1])
colnames(Amat) <- c(1:dim(Amat)[2])

# RHS constraints
zero_vector_RHS = rep(0,15)
RHS_model1 = c(1, 0.93, 0.07, 0.511, 0.106, 0.07, 0.03, 0)
RHS_model2 = c(1, 0.72, 0.07, 0.4355, 0.095, 0.05, 0.04, 0)
RHS_model3 = c(1, 0.53, 0.12, 0.2855, 0.063, 0.04, 0.06, 0.04)

# para_Zenic = c(6)
para_Yeats = c(8)
#RHS_model1_Zanic = c(RHS_model2, target2, zero_vector_RHS, para_Zenic)
RHS_model1_Yeats = c(RHS_model1, target1, zero_vector_RHS, para_Yeats)

# objective function 
obj1 = rep(0, 15)
obj2 = rep(3, 10)
obj3 = rep(2, 30)

Objvec = c(obj1, obj2, obj3, obj1)

#==============================================LP regression==============================================
## solve using glpk via the API
library(glpkAPI)

## initialize model
lp<- initProbGLPK()

## model dimensions
nrows <- dim(Amat)[1]
ncols <- dim(Amat)[2]

## use sparse format for model data
A_Zenic <- as.data.frame(as.table(Amat))
colnames(A_Zenic) <- c("row", "col", "value")

library(dplyr)
A_Zenic <- A_Zenic %>%filter(value != 0)

# load constraint matrix coefficients
nnonzero = dim(A_Zenic)[1]
rindexnonzero = A_Zenic$row
cindexnonzero = A_Zenic$col
valuesnonzero = A_Zenic$value

# row upper and lower bounds for 40 rows
lower1 = RHS_model1_Yeats[1:23]
rlower = c(lower1, rep(-1, 15), 7)
#rlower <- rep(-1000,nrows)
rupper <- RHS_model1_Yeats

# column upper and lower bounds for 70 columns
clower <- c(rep(0, 55),0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0)
cupper <- rep(1,ncols)

# maximize objective GLP_Max (minimize with GLP_MIN)
setObjDirGLPK(lp,GLP_MIN)

# tell model how many rows and columns
addRowsGLPK(lp, nrows)
addColsGLPK(lp, ncols)

# add column limits
setColsBndsGLPK(lp,c(1:ncols), clower, cupper)

# indicate column types
setColsKindGLPK(lp,1:ncols, c(rep(GLP_CV, 55), rep(GLP_BV, 15)))

# set row bounds
setRowsBndsGLPK(lp,c(1:nrows),rlower,rupper)

# set objective coefficients
setObjCoefsGLPK(lp,c(1:ncols),Objvec)

loadMatrixGLPK(lp,nnonzero, rindexnonzero, cindexnonzero,valuesnonzero)
# solve MIP problem and solve LP relaxation via simplex method first
solveSimplexGLPK(lp)

# add column limits
setColsBndsGLPK(lp,c(1:ncols), clower, cupper)

# indicate column types
setColsKindGLPK(lp,1:ncols, c(rep(GLP_CV, 55), rep(GLP_IV, 15)))

# set row bounds
setRowsBndsGLPK(lp,c(1:nrows),rlower,rupper)

solveMIPGLPK(lp)

# get retsults of MIP solution
mipStatusGLPK(lp)
status_codeGLPK(mipStatusGLPK(lp))

# report results
mipObjValGLPK(lp)
mipColsValGLPK(lp)
mipRowsValGLPK(lp)

getRowDualGLPK(lp,1)
getRowsDualGLPK(lp)