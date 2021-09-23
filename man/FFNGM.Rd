\name{FFNGM}
\alias{FFNGM}
\title{Fitting and Forecasting using Nonlinear Growth Models}
\usage{FFNGM (x, t, model=c("Monomolecular", "Logistic", "Gompertz"), k, y, r, h)}
\description{The FFNGM function fits nonlinear growth models to time series data and computes the h step ahead forecast values.
}

%-
\arguments{
  \item{x}{
a univariate time series data.
}
\item{t}{
a numeric vector containing time points.
}
\item{model}{
“Monomolecular” or “Logistic” or “Gompertz”.

}
\item{k}{
Initial estimate of carrying capacity (maximum limit of the considered time series data).

}
\item{y}{
Initial estimate of starting value of the considered time series data.

}
\item{r}{
Initial estimate of growth rate.}
\item{h}{
The forecast horizon.
}
}

\details{Using the nonlinear least squares method, this function estimates the parameters of nonlinear growth models for time series data. This function returns the fitted model summary, as well as the model's fitted values and various evaluation criteria. This function also returns the fitted model's h step ahead forecasted values.

}
\value{
\item{modelsummary }{Summary of the fitted model}
  \item{fitted.values }{Fitted  values of the model}
  \item{MAE}{Mean Absolute Error (MAE) of the fitted  model}
  \item{MAPE}{Mean Absolute Percentage Error (MAPE) of the fitted  model}
  \item{MSE}{Mean Square Error (MSE) of fitted the  model}
  \item{RMSE}{Root Mean Square Error (RMSE) of the fitted  model}
  \item{forecasted.values}{h step ahead forecasted values of the fitted  Model}
  }
\references{
%%
Pal, S. and Mazumder, D. (2015). Forecasting groundnut production of India using nonlinear growth models. Journal Crop and Weed, 11, 67-70.

Seber, G. A.F. and Wild, C. J. 2003. Nonlinear Regression, 2, New York: John Wiley.


}
\author{
Mrinmoy Ray, K. N. Singh, Kanchan Sinha, Rajeev Ranjan Kumar, Prakash Kumar
}

\seealso{nls}
\examples{
t<-c(0:12)
x<-c(57.97,66.02, 72.62, 77.87, 81.95, 85.07, 87.43, 89.20, 90.52, 91.50, 92.22, 92.75, 93.1)
FFNGM(x,t,"Gompertz",94, 55, 0.1, 3)

}

\keyword{nonlinear growth models}
\keyword{Time Series}
\keyword{Gompertz}
\keyword{Logistic}
\keyword{Monomolecular}

