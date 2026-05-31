# ==============================================================================
# parameter_estimates_Weibull.R
# МП-оценки параметров распределения Вейбулла
# Полная, правоцензурированная и интервально-цензурированная выборки
#
# Параметризация: F(t; beta, eta) = 1 - exp( -(t/eta)^beta )
#   beta — параметр формы
#   eta  — параметр масштаба
# ==============================================================================

# ------------------------------------------------------------------------------
# МПО по полной выборке (Ньютон–Рафсон)
#
# beta находится из уравнения правдоподобия:
#   h(beta) = 1/beta + mean(ln t) - sum(t^beta ln t) / sum(t^beta) = 0
# eta вычисляется аналитически:
#   eta_hat = (sum(t^beta_hat) / n)^{1/beta_hat}
#
# Аргументы:
#   data     — числовой вектор времён наблюдений (> 0)
#   beta0    — начальное приближение для beta
#   tol      — точность сходимости
#   max_iter — максимальное число итераций
# Возвращает: список (eta_hat, beta_hat, iterations, converged)
# ------------------------------------------------------------------------------

mle_weibull_complete <- function(data, beta0 = 1.0,
                                  tol = 1e-9, max_iter = 1000L) {
  data <- as.numeric(data)
  n    <- length(data)

  if (n == 0L) stop("Выборка пуста")
  if (any(data <= 0)) stop("Времена наблюдений должны быть положительными")

  ln_t      <- log(data)
  mean_ln_t <- mean(ln_t)

  h <- function(beta) {
    t_beta    <- data^beta
    t_beta_ln <- t_beta * ln_t
    1.0 / beta + mean_ln_t - sum(t_beta_ln) / sum(t_beta)
  }

  h_prime <- function(beta) {
    t_beta     <- data^beta
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
  eta_hat  <- (sum(data^beta_hat) / n)^(1.0 / beta_hat)

  list(
    eta_hat    = eta_hat,
    beta_hat   = beta_hat,
    iterations = iteration,
    converged  = converged
  )
}

# ------------------------------------------------------------------------------
# МПО по правоцензурированной выборке (I или II тип, Ньютон–Рафсон)
#
# Суммы t^beta берутся по ВСЕМ N наблюдениям,
# среднее ln t — только по n отказам.
#
# Аргументы:
#   times      — числовой вектор времён наблюдений t_1,...,t_N (> 0)
#   indicators — целочисленный вектор: 1 — отказ, 0 — цензурированное наблюдение
#   beta0      — начальное приближение для beta
#   tol        — точность сходимости
#   max_iter   — максимальное число итераций
# Возвращает: список (eta_hat, beta_hat, iterations, converged)
# ------------------------------------------------------------------------------

mle_weibull_censored <- function(times, indicators, beta0 = 1.0,
                                  tol = 1e-9, max_iter = 1000L) {
  times      <- as.numeric(times)
  indicators <- as.integer(indicators)
  N <- length(times)
  n <- sum(indicators)

  if (n == 0L) stop("Нет ни одного отказа — оценка невозможна")
  if (any(times <= 0)) stop("Времена наблюдений должны быть положительными")

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
# МПО по интервально-цензурированной выборке
# Двумерный Ньютон–Рафсон с численным Гессианом
#
# Аргументы:
#   t1       — левые границы интервалов (правая цензура: t2 = Inf; левая: t1 = 0)
#   t2       — правые границы интервалов
#   eta0     — начальное приближение для eta (NULL — автоматически)
#   beta0    — начальное приближение для beta (NULL — beta0 = 1)
#   tol      — точность сходимости
#   max_iter — максимальное число итераций
# Возвращает: список (eta_hat, beta_hat, iterations, converged, score_eta, score_beta)
# ------------------------------------------------------------------------------

mle_weibull_interval <- function(t1, t2, eta0 = NULL, beta0 = NULL,
                                  tol = 1e-9, max_iter = 1000L) {
  t1 <- as.numeric(t1)
  t2 <- as.numeric(t2)
  N  <- length(t1)

  if (N == 0L) stop("Выборка пуста")
  if (length(t2) != N) stop("Длины t1 и t2 должны совпадать")
  if (any(t1 < 0)) stop("Левые границы t1 должны быть неотрицательными")
  if (any(t2 <= t1)) stop("Должно выполняться t2 > t1")

  is_right    <- is.infinite(t2)
  is_left     <- (t1 == 0) & !is_right
  is_interval <- !is_right & !is_left

  if (is.null(eta0)) {
    midpoints <- ifelse(is_right, t1 * 2,
                        ifelse(is_left, t2 / 2, (t1 + t2) / 2))
    eta0 <- mean(midpoints)
    if (eta0 <= 0) eta0 <- 1.0
  }
  if (is.null(beta0)) beta0 <- 1.0

  compute_score <- function(eta, beta) {
    S_eta  <- 0.0
    S_beta <- 0.0

    if (any(is_right)) {
      t1_r      <- t1[is_right]
      a         <- (t1_r / eta)^beta
      S_eta     <- S_eta + sum(a)
      S_beta    <- S_beta + sum(a * log(t1_r / eta))
    }

    if (any(is_left)) {
      t2_l      <- t2[is_left]
      b         <- (t2_l / eta)^beta
      eb        <- exp(-b)
      denom     <- pmax(1.0 - eb, 1e-300)
      S_eta     <- S_eta + sum(-b * eb / denom)
      S_beta    <- S_beta + sum(b * log(t2_l / eta) * eb / denom)
    }

    if (any(is_interval)) {
      t1_iv     <- t1[is_interval]
      t2_iv     <- t2[is_interval]
      a         <- (t1_iv / eta)^beta
      b         <- (t2_iv / eta)^beta
      ea        <- exp(-a)
      eb        <- exp(-b)
      D         <- pmax(ea - eb, 1e-300)
      S_eta     <- S_eta + sum((a * ea - b * eb) / D)
      S_beta    <- S_beta + sum(
        (b * log(t2_iv / eta) * eb - a * log(t1_iv / eta) * ea) / D
      )
    }

    c(S_eta, S_beta)
  }

  # Численный Гессиан (центральные разности)
  compute_hessian <- function(eta, beta, h_rel = 1e-5) {
    H      <- matrix(0, 2, 2)
    h_eta  <- h_rel * eta
    h_beta <- h_rel * beta

    S_plus  <- compute_score(eta + h_eta, beta)
    S_minus <- compute_score(eta - h_eta, beta)
    H[, 1]  <- (S_plus - S_minus) / (2 * h_eta)

    S_plus  <- compute_score(eta, beta + h_beta)
    S_minus <- compute_score(eta, beta - h_beta)
    H[, 2]  <- (S_plus - S_minus) / (2 * h_beta)

    H
  }

  eta       <- eta0
  beta      <- beta0
  converged <- FALSE

  for (iteration in seq_len(max_iter)) {
    S <- compute_score(eta, beta)
    H <- compute_hessian(eta, beta)

    delta <- tryCatch(solve(H, S), error = function(e) S * 0.1)

    eta_new  <- eta  - delta[1]
    beta_new <- beta - delta[2]

    if (eta_new  <= 0) eta_new  <- eta  / 2.0
    if (beta_new <= 0) beta_new <- beta / 2.0

    if (abs(eta_new - eta) < tol && abs(beta_new - beta) < tol) {
      converged <- TRUE
      eta       <- eta_new
      beta      <- beta_new
      break
    }

    eta  <- eta_new
    beta <- beta_new
  }

  S_final <- compute_score(eta, beta)

  list(
    eta_hat    = eta,
    beta_hat   = beta,
    iterations = iteration,
    converged  = converged,
    score_eta  = S_final[1],
    score_beta = S_final[2]
  )
}
