# ==============================================================================
# right_censoring_Weibull.R
# Бутстрэп-критерий Жу для правоцензурированных данных II типа:
# проверка распределения Вейбулла
# Zhu (2020), Algorithm III
# ==============================================================================

# ------------------------------------------------------------------------------
# МПО параметров Вейбулла по правоцензурированным данным (Ньютон–Рафсон)
#
# Параметризация: F(t; beta, eta) = 1 - exp( -(t/eta)^beta )
#   beta — параметр формы
#   eta  — параметр масштаба
# ------------------------------------------------------------------------------

mle_weibull_censored <- function(times, indicators, beta0 = 1.0,
                                  tol = 1e-9, max_iter = 1000L) {
  times      <- as.numeric(times)
  indicators <- as.integer(indicators)
  N <- length(times)
  n <- sum(indicators)   # число отказов

  if (n == 0L) stop("Нет ни одного отказа — оценка невозможна")
  if (any(times <= 0)) stop("Времена должны быть положительными")

  ln_t               <- log(times)
  mean_ln_t_failures <- sum(indicators * ln_t) / n

  h <- function(beta) {
    t_beta    <- times^beta
    t_beta_ln <- t_beta * ln_t
    S0 <- sum(t_beta)
    S1 <- sum(t_beta_ln)
    1.0 / beta + mean_ln_t_failures - S1 / S0
  }

  h_prime <- function(beta) {
    t_beta     <- times^beta
    t_beta_ln  <- t_beta * ln_t
    t_beta_ln2 <- t_beta * ln_t^2
    S0 <- sum(t_beta)
    S1 <- sum(t_beta_ln)
    S2 <- sum(t_beta_ln2)
    -1.0 / beta^2 - (S2 * S0 - S1^2) / S0^2
  }

  beta      <- beta0
  converged <- FALSE

  for (iteration in seq_len(max_iter)) {
    h_val  <- h(beta)
    hp_val <- h_prime(beta)

    if (abs(hp_val) < 1e-15) break

    beta_new <- beta - h_val / hp_val
    if (beta_new <= 0) beta_new <- beta / 2.0

    if (abs(beta_new - beta) < tol) {
      converged <- TRUE
      beta      <- beta_new
      break
    }
    beta <- beta_new
  }

  beta_hat <- beta
  eta_hat  <- (sum(times^beta_hat) / n)^(1.0 / beta_hat)

  list(
    eta_hat    = eta_hat,
    beta_hat   = beta_hat,
    iterations = iteration,
    converged  = converged
  )
}

# ------------------------------------------------------------------------------
# Обёртка МПО для выборки II типа цензуры
# x_obs — наблюдаемые значения (первые m из n)
# ------------------------------------------------------------------------------

mle_weibull_typeII <- function(x_obs, n) {
  m          <- length(x_obs)
  times      <- c(x_obs, rep(x_obs[m], n - m))
  indicators <- c(rep(1L, m), rep(0L, n - m))
  mle_weibull_censored(times, indicators)
}

# ------------------------------------------------------------------------------
# Статистики Крамера–фон Мизеса и Андерсона–Дарлинга
# для выборки II типа цензуры (Pakyari & Balakrishnan, 2012)
# C_{i:m:n} = i/n при схеме R = (0,...,0, n-m)
# ------------------------------------------------------------------------------

stat_W2_typeII <- function(u, m, n) {
  mid <- (2 * seq_len(m) - 1) / (2 * n)
  sum((u - mid)^2) / n
}

stat_A2_typeII <- function(u, m, n) {
  u   <- pmin(pmax(u, 1e-15), 1 - 1e-15)
  mid <- (2 * seq_len(m) - 1) / (2 * n)
  sum((u - mid)^2 / (u * (1 - u))) / n
}

# ------------------------------------------------------------------------------
# Основная функция: бутстрэп-критерий Жу (2020), Algorithm III
#
#   H0: данные из распределения Вейбулла (параметры неизвестны)
#
# Аргументы:
#   x_obs  — наблюдаемые значения (отсортированные первые m из n)
#   n      — полный объём выборки
#   B      — число бутстрэп-реплик
#   alpha  — уровень значимости
# ------------------------------------------------------------------------------

zhu_weibull_test <- function(x_obs, n, B = 499, alpha = 0.05) {

  m     <- length(x_obs)
  x_obs <- sort(x_obs)

  # Шаг 2: МПО и статистики по наблюдённым данным
  fit0 <- mle_weibull_typeII(x_obs, n)
  if (!fit0$converged) stop("МПО не сошлось для исходных данных")

  b0 <- fit0$beta_hat
  e0 <- fit0$eta_hat

  u0   <- pweibull(x_obs, shape = b0, scale = e0)
  W2_0 <- stat_W2_typeII(u0, m, n)
  A2_0 <- stat_A2_typeII(u0, m, n)

  # Шаг 3: параметрический бутстрэп
  W2_boot <- numeric(B)
  A2_boot <- numeric(B)
  ok      <- logical(B)

  for (j in seq_len(B)) {
    x_star <- sort(rweibull(n, shape = b0, scale = e0))[seq_len(m)]

    fit_j <- tryCatch(
      mle_weibull_typeII(x_star, n),
      error = function(e) NULL
    )

    if (is.null(fit_j) || !fit_j$converged) next

    u_j <- pweibull(x_star, shape = fit_j$beta_hat, scale = fit_j$eta_hat)

    W2_boot[j] <- stat_W2_typeII(u_j, m, n)
    A2_boot[j] <- stat_A2_typeII(u_j, m, n)
    ok[j]      <- TRUE
  }

  W2_boot <- W2_boot[ok]
  A2_boot <- A2_boot[ok]
  B_ok    <- sum(ok)

  if (B_ok < 10) stop(sprintf("Слишком мало валидных бутстрэп-реплик: %d", B_ok))

  # Шаг 4: p-значения
  p_W2 <- (sum(W2_boot >= W2_0) + 1) / (B_ok + 1)
  p_A2 <- (sum(A2_boot >= A2_0) + 1) / (B_ok + 1)

  # Шаг 5: решение
  W2_reject <- p_W2 <= alpha
  A2_reject <- p_A2 <= alpha

  result <- list(
    test           = "Крамера–фон Мизеса и Андерсона–Дарлинга",
    distribution   = "Вейбулла",
    data_type      = "правоцензурированная выборка II типа",
    method         = "параметрический бутстрэп",
    hypothesis     = "H0: выборка согласуется с распределением Вейбулла",
    alpha          = alpha,
    eta_hat        = e0,
    beta_hat       = b0,
    n              = n,
    m              = m,
    censoring_rate = (n - m) / n,
    B              = B,
    B_valid        = B_ok,
    W2             = W2_0,
    A2             = A2_0,
    p_W2           = p_W2,
    p_A2           = p_A2,
    W2_reject      = W2_reject,
    A2_reject      = A2_reject
  )

  cat("Критерии согласия:", result$test, "\n")
  cat("Распределение:", result$distribution, "\n")
  cat("Тип данных:", result$data_type, "\n")
  cat("Метод:", result$method, "\n")
  cat(result$hypothesis, "\n\n")

  cat("Оценки параметров:\n")
  cat("  eta_hat  =", round(result$eta_hat, 4), "\n")
  cat("  beta_hat =", round(result$beta_hat, 4), "\n")
  cat("  alpha    =", round(result$alpha, 4), "\n\n")

  cat("Выборка:\n")
  cat("  n =", result$n, "\n")
  cat("  m =", result$m, "(наблюдаемых)\n")
  cat("  доля цензурирования =", round(result$censoring_rate, 4), "\n")
  cat("  бутстрэп-реплик (валидных) =", result$B_valid, "из", result$B, "\n\n")

  cat("Результаты:\n")
  cat("  Крамер–фон Мизес:\n")
  cat("    W2      =", round(result$W2, 4), "\n")
  cat("    p-value =", round(result$p_W2, 4), "\n")
  cat("    Вывод  ->",
      ifelse(result$W2_reject, "H0 отклоняется", "H0 не отклоняется"), "\n\n")

  cat("  Андерсон–Дарлинг:\n")
  cat("    A2      =", round(result$A2, 4), "\n")
  cat("    p-value =", round(result$p_A2, 4), "\n")
  cat("    Вывод  ->",
      ifelse(result$A2_reject, "H0 отклоняется", "H0 не отклоняется"), "\n\n")

  return(invisible(result))
}

# ==============================================================================
# Пример использования (не выполняется при source)
# ==============================================================================

if (FALSE) {

  set.seed(123)
  n     <- 100
  x_obs <- sort(rweibull(n, shape = 2, scale = 3))[1:80]  # m = 80 из n = 100
  zhu_weibull_test(x_obs, n = n, B = 499, alpha = 0.05)
}
