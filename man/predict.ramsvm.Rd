\name{predict}
\alias{predict}
\alias{predict,ramsvm-method}
\title{A function that provides class label predictions for objects 
       returned by the "ramsvm" function.}
\description{This function provides predictions on a test data set 
             using the obtained classifier from a call of the ramsvm 
             function.}
\usage{
\S4method{predict}{ramsvm}(object,newdata=NULL,lambda=NULL,...)
}
\arguments{
  \item{object}{An object returned by the "ramsvm" function.}
  \item{newdata}{The new predictor matrix. The number and order of 
                 predictors in newdata should be the same as those of 
                 x, which is used in the ramsvm function. If not 
                 specified, the program uses the training x matrix 
                 as the prediction object.}
  \item{lambda}{The lambda values on which the user wants to predict. 
                If not specified, the program will use the lambda 
                values returned by the ramsvm function.}
  \item{...}{Not used.}
 }

\value{
  \item{pred.y}{A prediction is made for each lambda  
                value.}    
}
\references{
C. Zhang, Y. Liu, J. Wang and H. Zhu. (2015+). 
Reinforced Angle-based Multicategory Support Vector Machines.
Journal of Computational and Graphical Statistics, in press.}
\author{Chong Zhang, Yufeng Liu, and Shannon Holloway}

\seealso{\code{\link{ramsvm}}} 

\examples{
data(iris)
a <- ramsvm(x = as.matrix(iris[,-5]),
            y = iris[,5],
            lambda = 0.2,
            kernel="gaussian")
predict(a)
}
