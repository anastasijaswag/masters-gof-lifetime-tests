# ==============================================================================
# parameter_estimates_Exp.R
# МП-оценки параметра масштаба экспоненциального распределения
# Полная, правоцензурированная и интервально-цензурированная выборки
# ==============================================================================

# ------------------------------------------------------------------------------
# МПО по полной выборке
#
# F(t) = 1 - exp(-t/eta),  eta_hat = mean(t)
#
# Аргументы:
#   data — числовой вектор времён наблюдений t_1, ..., t_n
# Возвращает: eta_hat (скаляр)
# ------------------------------------------------------------------------------

mle_exp_complete <- function(data) {
  data <- as.numeric(data)
  n    <- length(data)

  if (n == 0L) stop("Выборка пуста")
  if (any(data < 0)) stop("Времена наблюдений должны быть неотрицательными")

  mean(data)
}

# ------------------------------------------------------------------------------
# МПО по правоцензурированной выборке (I или II тип)
#
# eta_hat = sum(t_i) / n_failures
#
# Аргументы:
#   times      — числовой вектор времён наблюдений (отказы и моменты цензуры)
#   indicators — целочисленный вектор: 1 — отказ, 0 — цензурированное наблюдение
# Возвращает: eta_hat (скаляр)
# ------------------------------------------------------------------------------

mle_exp_censored <- function(times, indicators) {
  times      <- as.numeric(times)
  indicators <- as.integer(indicators)
  N <- length(times)
  n <- sum(indicators)   # число отказов

  if (N == 0L) stop("Выборка пуста")
  if (n == 0L) stop("Нет ни одного отказа — оценка невозможна")
  if (any(times < 0)) stop("Времена наблюдений должны быть неотрицательными")

  sum(times) / n
}

# ------------------------------------------------------------------------------
# МПО по интервально-цензурированной выборке (Ньютон–Рафсон)
#
# Скор-функция: S(eta) = sum_i g_i(eta) = 0
# g_i(eta) зависит от типа наблюдения:
#   правая цензура: g_i = t_{1i}
#   левая цензура:  g_j = -t_{2j} exp(-t_{2j}/eta) / (1 - exp(-t_{2j}/eta))
#   интервал:       g_i = (t_{1i} e1 - t_{2i} e2) / (e1 - e2)
#
# Аргументы:
#   t1       — левые границы интервалов (правая цензура: t2 = Inf; левая: t1 = 0)
#   t2       — правые границы интервалов
#   eta0     — начальное приближение (NULL — автоматически)
#   tol      — точность сходимости
#   max_iter — максимальное число итераций
# Возвращает: список (eta_hat, iterations, converged, score)
# ------------------------------------------------------------------------------

mle_exp_interval <- function(t1, t2, eta0 = NULL,
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
    midpoints <- ifelse(is_right, t1 * 2, (t1 + t2) / 2)
    eta0 <- mean(midpoints)
    if (eta0 <= 0) eta0 <- 1.0
  }

  compute_terms <- function(eta) {
    g <- numeric(N)

    if (any(is_right)) {
      g[is_right] <- t1[is_right]
    }

    if (any(is_left)) {
      t2_l <- t2[is_left]
      e2   <- exp(-t2_l / eta)
      g[is_left] <- -t2_l * e2 / (1.0 - e2)
    }

    if (any(is_interval)) {
      t1_iv <- t1[is_interval]
      t2_iv <- t2[is_interval]
      e1 <- exp(-t1_iv / eta)
      e2 <- exp(-t2_iv / eta)
      g[is_interval] <- (t1_iv * e1 - t2_iv * e2) / (e1 - e2)
    }

    g
  }

  score_fun <- function(eta) sum(compute_terms(eta))

  score_deriv <- function(eta) {
    g  <- compute_terms(eta)
    dg <- numeric(N)

    if (any(is_left)) {
      t2_l   <- t2[is_left]
      e2     <- exp(-t2_l / eta)
      ratio  <- -t2_l^2 * e2 / (1.0 - e2)
      dg[is_left] <- (1.0 / eta^2) * (ratio - g[is_left]^2)
    }

    if (any(is_interval)) {
      t1_iv <- t1[is_interval]
      t2_iv <- t2[is_interval]
      e1    <- exp(-t1_iv / eta)
      e2    <- exp(-t2_iv / eta)
      ratio <- (t1_iv^2 * e1 - t2_iv^2 * e2) / (e1 - e2)
      dg[is_interval] <- (1.0 / eta^2) * (ratio - g[is_interval]^2)
    }

    sum(dg)
  }

  eta       <- eta0
  converged <- FALSE

  for (iteration in seq_len(max_iter)) {
    S_val <- score_fun(eta)
    S_der <- score_deriv(eta)

    if (abs(S_der) < 1e-15) break

    eta_new <- eta - S_val / S_der
    if (eta_new <= 0) eta_new <- eta / 2.0

    if (abs(eta_new - eta) < tol) {
      converged <- TRUE
      eta       <- eta_new
      break
    }
    eta <- eta_new
  }

  list(
    eta_hat    = eta,
    iterations = iteration,
    converged  = converged,
    score      = score_fun(eta)
  )
}
