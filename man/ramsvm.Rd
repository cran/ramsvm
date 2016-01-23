\name{ramsvm}
\alias{ramsvm}
\title{The classifier for Reinforced Angle-Based Multicategory 
       Support Vector Machines (RAMSVMs).}
\description{A function that provides the RAMSVMs classifier for 
             linear learning, polynomial learning, and kernel learning.}
\usage{
ramsvm(x, y, lambda, gamma = 0.5, weight = NULL, kernel = "linear", 
       kparam = NULL, large = FALSE, epsilon = NULL, warm = NULL, 
       nb.core=NULL)
}
\arguments{
  \item{x      }{The x matrix for the training dataset. 
                 Columns represent the covariates, and rows represent 
                 the instances. There should be no NA/NaN values in x.}
  \item{y      }{The labels for the training dataset.}
  \item{gamma  }{The convex combination parameter of the loss function.}
  \item{weight }{The weight vector for each observation. By default, 
                 the program uses equal weights for all observations.}
  \item{lambda }{The user specified lambda values.}
  \item{kernel }{The kernel for classification.}
  \item{kparam }{The kernel parameter. If kernel=linear, this option 
                 is ignored. For kernel=polynomial, it is the order 
                 of the polynomial functions. For kernel=gaussian, 
                 it is the Gaussian kernel parameter.}
  \item{large  }{Whether the number of observations is large in the 
                 data. If TRUE, then the algorithm will split the 
                 data set into several parts and train on each part 
                 to provide a warm start for the entire data training. 
                 This option aims to enhance the computational speed.}
  \item{epsilon}{Convergence threshold in coordinate descent circling 
                 algorithm. The smaller epsilon is, the more accurate 
                 the final model is, and the more time it takes for 
                 calculation. Default is 
                 (0.0001*number of observations*number of classes).}
  \item{warm   }{A matrix that contains the warm start for slack 
                 variables alpha. This option is especially useful 
                 when the user wishes to obtain the classifier with 
                 higher level accuracy (smaller epsilon) or with a 
                 different lambda, if the warm start is available 
                 from an existing ramsvm output.}
  \item{nb.core}{The number of threads to use for parallel computing. 
                 If null, the code will automatically detect and use 
                 the number of CPU cores. This option is used only 
                 when large=TRUE.}
}

\value{
  An object of class ramsvm is returned.

  If kernel=linear, this S4 object contains the following:
  \item{x      }{A copy of the input covariate matrix.}
  \item{y      }{A copy of the input labels.}
  \item{y.name }{The class names of y.}
  \item{k      }{Number of classes in the classification problems.}
  \item{gamma  }{A copy of the convex combination parameter of the loss 
                 function.}
  \item{weight }{The weight vector for each observation.}
  \item{lambda }{The lambda vector of all lambdas in the solution path.}
  \item{beta   }{A list of matrices containing the estimated parameters 
                 of the classification function. Each matrix in the 
                 list corresponds to the lambda value in the solution 
                 path in order. For one single matrix, the rows 
                 correspond to a specific predictor.}
  \item{beta0  }{A list of the intercepts of the classification 
                 function. Each vector in the list corresponds to 
                 the lambda in the solution path in order.}
  \item{epsilon}{Convergence threshold in coordinate descent circling 
                 algorithm.}
  \item{call}{The call of ramsvm.}

  If kernel != linear the S4 object also contains the following:

  \item{kernel }{The kernel for classification.}
  \item{kparam }{The kernel parameter.}
}

\references{
C. Zhang, Y. Liu, J. Wang and H. Zhu. (2015+). 
Reinforced Angle-based Multicategory Support Vector Machines.
Journal of Computational and Graphical Statistics, in press.}
\author{Chong Zhang, Yufeng Liu, and Shannon Holloway}

\seealso{\code{\link{predict}}} 

\examples{
data(iris)
ramsvm(x = as.matrix(iris[,-5]),
       y = iris[,5],
       lambda = 0.2,
       kernel="gaussian")
}
