#' Parameter One-Sample T-Test
#' @param data data yang digunakan untuk uji one sample t test
#' @param mu_0 rata-rata populasi yang dijadikan acuan dalam pengujian hipotesis
#' @return nilai t hitung
#' @examples
#' data <- c(140,155,160,145,150,138,152,148,157,143)
#' mu_0 <- 150
#' @export

#Uji Satu Arah (One-Sample T-Test)
one_sample_t_test <- function(data, mu_0) {
  x_bar <- mean(data)
  s <- sd(data)
  n <- length(data)
  result <- (x_bar - mu_0) / (s / sqrt(n))
  return(result)
}

#' Parameter Two-Sample T-Test
#' @param data1 data pertama yang digunakan untuk uji Two sample t test
#' @param data2 data kedua yang digunakan untuk uji Two sample t test
#' @return nilai t hitung
#' @examples
#' data1 <- c(78,85,88,92,76,81,90,87,84,79)
#' data2 <- c(72,80,75,78,74,79,77,81,76,73)
#' @export

#Uji Dua Arah (Two-Sample T-Test)
two_sample_t_test <- function(data1, data2) {
  x_bar1 <- mean(data1)
  x_bar2 <- mean(data2)
  s1_sq <- var(data1)
  s2_sq <- var(data2)
  n1 <- length(data1)
  n2 <- length(data2)
  result <- (x_bar1 - x_bar2) / sqrt((s1_sq / n1) + (s2_sq / n2))
  return(result)
}

#' Parameter Uji Linier Regresi
#' @param x Vektor numerik yang berisi nilai variabel independen (prediktor)
#' @param y Vektor numerik yang berisi nilai variabel dependen (respon)
#' @return Intercept, Slope, R Squared, P-Value
#' example
#' x <- c(1, 2, 3, 4, 5)
#' y <- c(2, 4, 5, 4, 5)
#' @export

#Uji Linier
linear_regression_test <- function(x, y) {

x_bar <- mean(x)
y_bar <- mean(y)

# Menghitung Slope (beta_1)
atas <- sum((x - x_bar) * (y - y_bar))
bawah <- sum((x - x_bar)^2)
beta_1 <- atas / bawah

# Menghitung Intercept (beta_0)
beta_0 <- y_bar - beta_1 * x_bar

# Menghitung Koefisien Determinasi
y_hat <- beta_0 + beta_1 * x
residuals <- y - y_hat
ss_res <- sum(residuals^2)
ss_tot <- sum((y - y_bar)^2)
r_squared <- 1 - (ss_res / ss_tot)

# Menghitung P-value
se_beta_1 <- sqrt(ss_res / ((length(x) - 2) * sum((x - x_bar)^2)))
t_statistic <- beta_1 / se_beta_1
df <- length(x) - 2
p_value <- 2 * pt(-abs(t_statistic), df)
  result <- list (Intercept = beta_0, slope = beta_1, R_squared = r_squared, P_Value = p_value)
  return(result)
}


#' Parameter Uji regresi logistik
#' @param x Vektor numerik yang berisi nilai variabel independen (prediktor)
#' @param y Vektor biner (0/1) yang berisi nilai variabel dependen (respon)
#' @return Intercept, Slope, R Squared, P-Value
#' example
#' x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' y <- c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' @export

# Fungsi untuk regresi logistik
logistic_regression_test <- function(x, y) {

  # Menambahkan intercept ke variabel prediktor
  x <- cbind(1, x)  # Menambahkan kolom 1 untuk intercept
  n <- length(y)

  # Fungsi Log-Likelihood
  log_likelihood <- function(beta) {
    p_hat <- 1 / (1 + exp(-(x %*% beta)))  # Fungsi sigmoid untuk probabilitas
    ll <- sum(y * log(p_hat) + (1 - y) * log(1 - p_hat))
    return(-ll)  # Negatif karena kita ingin memaksimalkan likelihood
  }

  # Estimasi koefisien menggunakan optimasi MLE
  initial_beta <- rep(0, ncol(x))  # Inisialisasi beta dengan nol
  optim_result <- optim(initial_beta, log_likelihood)  # Maksimalkan log-likelihood

  # Estimasi koefisien beta (intercept dan slope)
  beta_hat <- optim_result$par

  # Menghitung prediksi (probabilitas)
  p_hat <- 1 / (1 + exp(-(x %*% beta_hat)))

  # Pseudo R-squared (McFadden)
  log_likelihood_model <- -optim_result$value
  log_likelihood_null <- sum(y * log(mean(y)) + (1 - y) * log(1 - mean(y)))
  pseudo_r_squared <- 1 - log_likelihood_model / log_likelihood_null

  # Menghitung standar error (menggunakan estimasi varians koefisien)
  hessian_matrix <- optim_result$hessian  # Matriks Hessian (turunan kedua dari log-likelihood)
  se_beta <- sqrt(diag(solve(-hessian_matrix)))  # Standar error dari beta

  # Menghitung Wald statistic dan p-value untuk beta_1
  t_statistic <- beta_hat[2] / se_beta[2]  # t-statistik untuk koefisien beta_1
  p_value <- 2 * (1 - pnorm(abs(t_statistic)))  # p-value dari distribusi normal

  # Hasil regresi logistik
  result <- list(Intercept = beta_hat[1], Slope = beta_hat[2],
                 Pseudo_R_squared = pseudo_r_squared, P_value = p_value)
  return(result)
}
