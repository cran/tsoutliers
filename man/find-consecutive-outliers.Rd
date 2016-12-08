\name{find.consecutive.outliers}
\alias{find.consecutive.outliers}

\title{Find outliers at consecutive time points}

\description{
Find outliers at consecutive time points and discard those with 
lower \eqn{t}-statistic in absolute value.
}

\usage{
find.consecutive.outliers(x, type)
}

\arguments{
\item{x}{a data frame containing the type, index and \eqn{t}-statistic related 
to each outlier, named respectively by columns: \code{type}, \code{ind} and \code{tstat}}.
\item{type}{a character, the type of outlier to be checked for runs at consecutive points.}
}

\details{
In the procedure, this function is applied by type of outliers, 
e.g., two or more consecutive LS are reduced to one LS.
It is still possible to get two or more consecutive outliers of different type 
at consecutive time points. 
Alternatively consecutive outliers of any type could be reduced to a single outlier,
For the time being, this is kept in this way so that for example, the following sequence 
AO1, LS2, LS3,LS4 can collapse to AO-1,LS-4, i.e., the AO is kept.
}

\value{
The row names of input \code{x} pointing to the outliers to be discarded, 
i.e., outliers with lower significance than adjacent outliers of the same type.
}

\note{
This function is intended for internal use; 
no check is done for correctness of the arguments passed as input.
}

\keyword{ts}
