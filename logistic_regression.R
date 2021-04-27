# Linear Discriminant Analysis
# (the multiclass version of logistic regression)
# Author: Abdulrahman Sawalma


# This method is more widely used when we have mroe than two response classes
# it is more stable, especially if the distribution of the predictors is normal
# Reading is from the book: https://faculty.marshall.usc.edu/gareth-james/ISL/ISLR%20Seventh%20Printing.pdf
# Also, check this for ROC curve: https://towardsdatascience.com/understanding-the-roc-and-auc-curves-a05b68550b69
# using Bayes theorim (see the book) we can compute the probability of Y=k | X=x
# using the prior probability (pi_k) and f(x). pi_k is the overall probability that a randomly
# chosen observation comes from kth class. f(x) is the probability of X=x|Y=y
# This prability (pi_k) and function (f(x)) will be used to estimate the
# posterior probabiltiy (probability for each observation)

# Based on the Bayes classifier, we assign an observation to the class for which
# the p_k(x) is highest in the discriminant function
# LDA uses estimates for the mean for each class, prior probability for each class
# and sd for all classes to approximate the discriminant function. All what remains
# is to plug in the x value into the equation 4.17
# The boundary beyond which a variable is considered to belong to one of two classes
# is given by the equation (mu_1+mu_2)/2 ... for multivariate version, things are different

# In the multivariate case, the equation is the same, but the sd and mean is replaced with
# a covariance matrix and vector of means, respectively.
# The quadratic Discriminant Analysis assumes that each of the classes has different 
# covariance matrix
# Shortly, if the classes shar common variance, LDA is better. If this assumption is far 
# off, QDA is better


# Let's start with the examples
library(ISLR)

# get to know the data
?Smarket
summary(Smarket)

#find the nrows, ncolumns of Smarket
dim(Smarket)

#pairwise correlations of all variables
cor(Smarket[,-9])

#produce a testing data
train = Smarket$Year<2005
Smarket.2005 = Smarket[!train,]
Direction.2005 = Smarket$Direction[!train]

#First, let's do logistic regression
glm.fits=glm(Direction ~ Lag1+Lag2+Lag3+Lag4+Lag5+Volume, data=Smarket,
             family=binomial, subset = train)
summary(glm.fits)

# Predict theprobabilities
glm.probs = predict(glm.fits, Smarket.2005, type = "response")
glm.probs[1:10]

#Notice that the given values are the probability of (Y=1 | X).
# We can know what does "1" mean by looking at the contrasts(Smarket$Direction)
# where we see that Up was given a value of 1, and Down of 0

# calculate accuracy of prediction
glm.pred=rep("Down",dim(Smarket.2005)[1])
glm.pred[glm.probs >.5]="Up"

# For inspection
table(glm.pred,Direction.2005)

print (paste0("Percentage of correct = ", 100 * mean(Direction.2005 == glm.pred),"%"))



# Now the LDA
library(MASS)
lda.fit = lda (Direction ~ Lag1+Lag2, data = Smarket, subset = train)
lda.fit

#######################################################################################
# Cross validation
#######################################################################################

library(ISLR)
library(boot)

set.seed(17)
glm.fit=glm(mpg~horsepower,data=Auto)
cv.glm(Auto ,glm.fit ,K=10)$delta[1]

Data = read.csv("/home/asawalma/git/data_analysis/cv_trial.csv")
Data$y = factor(Data$y)

equation = "Data$y ~ "
for (i in 1:length(colnames(Data[,-length(Data)]))){
  item = colnames(Data[,-length(Data)])[i]
  if (i==1){
    equation = paste0(equation,"Data$",item)
  }else{
    equation = paste0(equation," * Data$",item)
  }
}
equation = as.formula(equation)

glm.fit=glm(formula = equation, family = binomial())

cv.glm(data = Data ,glm.fit ,K=10)$delta[1]
