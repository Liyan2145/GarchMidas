#' This function estimates a multiplicative mixed-frequency GARCH model. For the sake of numerical stability, it is best to multiply log returns by 100.
#' @param data data frame containing a column named date of type 'Date'.
#' @param y name of high frequency dependent variable in df.
#' @param x covariate employed in GarchMidas.
#' @param K an integer specifying lag length K in the long-term component.
#' @param freq a string of the low frequency variable in the df.
#' @param GJR if TRUE, an asymmetric GJR-GARCH is used as the short-term component. If FALSE, a simple GARCH(1,1) is employed.
#' @keywords fit_GarchMidas
#' @export
#' @importFrom pracma jacobian
#' @importFrom stats nlminb
#' @importFrom pracma hessian
#' @importFrom pracma zeros
#' @importFrom stats constrOptim
#' @importFrom stats na.exclude
#' @importFrom stats optim
#' @importFrom stats pnorm
#' @importFrom stats var
#' @importFrom stats aggregate
#' @importFrom utils tail
#' @examples
#' \dontrun{
#' fit_GarchMidas(data = mu, y = "return", x = "epu", K = 24, freq = "month", GJR = TRUE)
#' }

fit_CarrMidas <- function(data, y, x, K, freq = "month", GJR = FALSE) {
  date_backup <- data[["date"]]
  data["date"] <- as.numeric(unlist(data["date"]))
  df <- data[, c(y, x, freq)]
  df[, freq] <- as.integer(unlist(df[ , freq]))
  g0 <- var(unlist(data[[y]]))
  covariate <- unlist(unique(data[c(freq, x)])[x])
  
  # Parameter estimation ----------------------------------------------------------------------------
  if(GJR){
    lf <- function(p) {
      mle(df = df, y = data[[y]], x = covariate, freq = freq, v = p["v"], 
          omega = 1 - p["alpha"] - p["beta"], alpha  = p["alpha"], beta = p["beta"], gamma = p["gamma"],
          m = p["m"], theta1 = p["theta1"], theta2 = p["theta2"], w1 = p["w1"], w2 = p["w2"], h0 = h0, K = K)
    }
    par0 <- c(v = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0, theta1 = 0, theta2 = 0, w1 = 1.00000001, w2 = 3)
    ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0, 0),
                    c(0,  0,  0, 0, 0, 0, 1, 0),
                    c(0,  0,  0, 0, 0, 0, 0, 1),
                    c(0,  1,  0, 0, 0, 0, 0, 0),
                    c(0,  0,  1, 0, 0, 0, 0, 0))
  }else{
    lf <- function(p) {
      mle(df = df, y = data[[y]], x = covariate, freq = freq, v = p["v"], 
          omega = 1 - p["alpha"] - p["beta"], alpha  = p["alpha"], beta = p["beta"], gamma = 0,
          m = p["m"], theta1 = p["theta1"], theta2 = p["theta2"], w1 = p["w1"], w2 = p["w2"], g0 = g0, K = K)
    }
    par0 <- c(v = 0, alpha = 0.02, beta = 0.85, m = 0, theta1 = 0, theta2 = 0, w1 = 1.00000001, w2 = 3)
    ui.opt <- rbind(c(0, -1, -1, 0, 0, 0, 0),
                    c(0,  0,  0, 0, 0, 1, 0),
                    c(0,  0,  0, 0, 0, 0, 1),
                    c(0,  1,  0, 0, 0, 0, 0),
                    c(0,  0,  1, 0, 0, 0, 0))
  }
  ci.opt <- c(-0.99999999, 1, 1, 0, 0)
  optim.par <- constrOptim(theta1 = par0, theta2 = par0, f = function(theta1, theta2) {theta <- c(theta1, theta2) sum(lf(theta))},
                            grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)
  par <- optim.par$par
  tau <- calculate_tau(df = data, x = covariate, freq = freq,
                          w1 = par["w1"], w2 = par["w2"],
                          theta1 = par["theta1"],theta2 = par["theta2"],
                          m = par["m"], K = K)$tau
  R2 <- unlist(data[y])
  if(GJR){
    gamma = par["gamma"]
  }else{
    gamma = 0
  }
  h <- c(rep(NA, times = sum(is.na((R2)/tau))),
         calculate_h(omega = 1 - par["alpha"] - par["beta"],
                     alpha = par["alpha"],
                     beta = par["beta"],
                     gamma = gamma,
                     as.numeric(na.exclude((R2)/tau)),
                     g0 = g0))
  df.fitted <- cbind(data[c("date", y, freq, x)], h = h, tau = tau)
  df.fitted$residuals <- unlist((df.fitted[y] - par["mu"]) / sqrt(df.fitted$g * df.fitted$tau))
  df.fitted$date <- as.Date(date_backup)
  
  # Standard errors --------------------------------------------------------------------------------
  # inv_hessian <- try({solve(-hessian(function (theta) {sum(lf(theta))},par,h=1e-6))})
  # if (class(inv_hessian) == "try-error") {
  #   warning("Inverting the Hessian matrix failed. Possible workaround: Multiply returns by 100.")
  #   inv_hessian<-zeros(length(par))
  # }
  # rob.std.err <- sqrt(diag(inv_hessian %*% crossprod(jacobian(lf, par)) %*% inv_hessian))
  # Output -----------------------------------------------------------------------------------------
  output <-
    list(parameters = par,
         # std.err = rob.std.err,
         residuals=df.fitted$residuals,
         # estimation = data.frame(estimate = round(par,3),
         #                         std.err = round(rob.std.err,3),
         #                         p.value = round(2 * (1 - pnorm(unlist(abs(par/rob.std.err)))),3)),
         phi=calculate_phi(w1 = par["w1"], w2 = par["w2"], K = K),
         tau = tau,
         h = h,
         df.fitted = df.fitted,
         K = K,
         llh = -optim.par$value,
         bic = log(sum(!is.na(tau))) * length(par) - 2 * (-optim.par$value),
         optim = optim.par)
  output$variance.ratio <- 100 *
                   var(log(aggregate(df.fitted$tau, by = df.fitted[freq],
                                     FUN = mean)[,2]),
                       na.rm = TRUE) /
                   var(log(aggregate(df.fitted$tau * df.fitted$h, by = df.fitted[freq],
                                     FUN = mean)[,2]),
                       na.rm = TRUE)
  class(output) <- "CarrMidas"
  output
}
