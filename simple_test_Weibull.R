# ==============================================================================
# simple_test_Weibull.R
# LB-GOF критерий: простая гипотеза H0: F = Weibull(beta, eta), оба параметра известны
# Ren (2003), Scand. J. Statist., Section 3, теоремы 1–2, ур. (7)–(13)
#
# Параметризация:
#   F(x; beta, eta) = 1 - exp( -(x/eta)^beta ),   x >= 0
#   beta — параметр формы
#   eta  — параметр масштаба
# ==============================================================================

source("utils_interval.R")

# ------------------------------------------------------------------------------
# Статистика Крамера–фон Мизеса для простой гипотезы (Вейбулл)
#
#   T*_m = 1/(12m) + sum_{i=1}^{m} [ U_{(i)} - (2i-1)/(2m) ]²
#
# где U_i = F0(X*_i) = pweibull(X*_i, beta, eta),  U_{(i)} — порядковые статистики
# (Shorack & Wellner, 1986, с. 147)
# ------------------------------------------------------------------------------

cvm_stat_simple_weibull <- function(x_star, beta, eta) {
  m <- length(x_star)
  u <- sort(pweibull(x_star, shape = beta, scale = eta))
  1 / (12 * m) + sum((u - (2 * seq_len(m) - 1) / (2 * m))^2)
}

# ------------------------------------------------------------------------------
# Выбор m (Ren, 2003, ур. (16)–(18))
# Специфика распределения входит только через F0, f0 и x_max
# ------------------------------------------------------------------------------

choose_m_simple_weibull <- function(fit, n, beta, eta,
                                    alpha     = 0.05,
                                    epsilon   = 0.02,
                                    eta_param = 0.10,
                                    gamma     = 1/3,
                                    K_pilot   = 30,
                                    delta_cap = 0.05) {

  # e = C_alpha - C_{alpha+epsilon}  (ур. (17))
  C_alpha     <- qCvM(1 - alpha)
  C_alpha_eps <- qCvM(1 - alpha - epsilon)
  e <- C_alpha - C_alpha_eps

  x_max  <- qweibull(1 - 1e-6, shape = beta, scale = eta)
  x_grid <- seq(0, x_max, length.out = 5000)

  F0_vals <- pweibull(x_grid, shape = beta, scale = eta)
  f0_vals <- dweibull(x_grid, shape = beta, scale = eta)
  Fn_vals <- stepfun_npmle(fit, x_grid)

  # r_n = int (F_n - F0)^2 dF0  (ур. (17))
  integrand_r <- (Fn_vals - F0_vals)^2 * f0_vals
  r_n <- trapez(x_grid, integrand_r)

  m0 <- max(2L, floor(n^gamma))

  # Оценка sigma_mn по пилотным LB-выборкам размера m = n^gamma
  cross_terms <- numeric(K_pilot)
  for (k in seq_len(K_pilot)) {
    x_star   <- sample_from_npmle(fit, m0)
    Fnm_vals <- ecdf(x_star)(x_grid)

    integrand_cross <-
      sqrt(m0) * (Fnm_vals - Fn_vals) *
      (n^gamma) * (Fn_vals - F0_vals) * f0_vals

    cross_terms[k] <- trapez(x_grid, integrand_cross)
  }

  sigma_mn <- sd(cross_terms)
  if (!is.finite(sigma_mn)) sigma_mn <- 0

  z_eta2 <- qnorm(1 - eta_param / 2)

  if (r_n < 1e-14) {
    m_hat <- floor(n^(2 * gamma - delta_cap))
  } else {
    val1  <- e / r_n
    denom <- sqrt(sigma_mn^2 * z_eta2^2 + e * r_n * n^(2 * gamma)) +
      sigma_mn * z_eta2
    val2  <- if (denom > 1e-14) {
      (e / denom)^2 * n^(2 * gamma)
    } else {
      n^(2 * gamma - delta_cap)
    }
    m_hat <- floor(min(val1, val2))
  }

  # m = o(n^{2gamma}), ограничиваем строго ниже n^{2gamma}
  cap <- max(2L, floor(n^(2 * gamma - delta_cap)))
  m   <- max(m0, m_hat)
  m   <- min(m, cap)
  m   <- max(2L, m)

  m
}

# ------------------------------------------------------------------------------
# Основная функция: LB-GOF критерий (простая гипотеза, Вейбулл)
#
#   H0 : F = Weibull(beta, eta)   (оба параметра известны)
#
# Аргументы:
#   L, R  — левые и правые границы интервалов наблюдения;
#            левая цензура: L = 0;  правая цензура: R = Inf
#   beta  — известный параметр формы (beta > 0)
#   eta   — известный параметр масштаба (eta > 0)
# ------------------------------------------------------------------------------

lb.gof.simple.weibull <- function(L, R, beta, eta,
                                  alpha     = 0.05,
                                  rho       = NULL,
                                  gamma     = 1/3,
                                  epsilon   = 0.02,
                                  eta_param = 0.10,
                                  K_pilot   = 30,
                                  n_pilot   = 10000,
                                  N_max     = 100000,
                                  N_min     = 30,
                                  verbose   = TRUE) {

  if (is.null(rho)) rho <- alpha / 2

  n <- length(L)
  stopifnot(length(R) == n, all(L <= R), beta > 0, eta > 0)

  # Шаг 1: NPMLE
  fit <- icfit(L, R)

  # Шаг 2: выбор m (ур. (16)–(18))
  m       <- choose_m_simple_weibull(
    fit = fit, n = n, beta = beta, eta = eta,
    alpha = alpha, epsilon = epsilon,
    eta_param = eta_param, gamma = gamma,
    K_pilot = K_pilot
  )
  C_alpha <- qCvM(1 - alpha)

  # Шаг 3: оценка p_n = P_n{ T*_m >= C_alpha }  (ур. (11))
  T_pilot <- numeric(n_pilot)
  for (k in seq_len(n_pilot)) {
    x_star     <- sample_from_npmle(fit, m)
    T_pilot[k] <- cvm_stat_simple_weibull(x_star, beta, eta)
  }
  p_n <- mean(T_pilot >= C_alpha)

  # Шаг 4: выбор N (ур. (12))
  N <- choose_N_lb(p_n = p_n, alpha = alpha, rho = rho,
                   N_max = N_max, N_min = N_min)$N

  # Шаг 5: внешний цикл LB-GOF (ур. (11), (13))
  W_j <- integer(N)
  for (j in seq_len(N)) {
    x_star  <- sample_from_npmle(fit, m)
    T_star  <- cvm_stat_simple_weibull(x_star, beta, eta)
    W_j[j]  <- as.integer(T_star >= C_alpha)
  }
  W_bar <- mean(W_j)

  # Шаг 6: правило отклонения (ур. (13))
  z_alpha_rho <- qnorm(1 - (alpha - rho))
  threshold   <- alpha + z_alpha_rho * sqrt(alpha * (1 - alpha) / N)
  reject      <- (W_bar >= threshold)

  if (verbose) {
    cat("Критерий согласия: Крамера–фон Мизеса (LB-GOF)\n")
    cat("Распределение: Вейбулла\n")
    cat("Тип данных: интервально-цензурированная выборка\n")
    cat("H0: F = Weibull(beta, eta), параметры известны\n\n")

    cat("Параметры гипотезы:\n")
    cat("  beta  =", round(beta, 4), "\n")
    cat("  eta   =", round(eta, 4), "\n")
    cat("  alpha =", round(alpha, 4), "\n\n")

    cat("Выборка:\n")
    cat("  n =", n, "\n\n")

    cat("Параметры метода:\n")
    cat("  m       =", m, "(размер LB-подвыборки)\n")
    cat("  N       =", N, "(число внешних повторений)\n")
    cat("  C_alpha =", round(C_alpha, 6), "\n\n")

    cat("Результат:\n")
    cat("  W_bar =", round(W_bar, 6), "\n")
    cat("  порог =", round(threshold, 6), "\n")
    cat("  Вывод ->",
        ifelse(reject, "H0 отклоняется", "H0 не отклоняется"), "\n\n")
  }

  list(
    reject    = reject,
    W_bar     = W_bar,
    threshold = threshold,
    m         = m,
    N         = N,
    p_n       = p_n,
    C_alpha   = C_alpha,
    alpha     = alpha,
    rho       = rho,
    beta      = beta,
    eta       = eta,
    npmle_fit = fit,
    mode      = "lb-gof-simple-weibull"
  )
}

# ==============================================================================
# Пример использования (не выполняется при source)
# ==============================================================================

if (FALSE) {

  set.seed(42)
  n    <- 200
  beta0 <- 2
  eta0  <- 3

  X <- rweibull(n, shape = beta0, scale = eta0)
  Y <- runif(n, 0, 6)
  Z <- Y + runif(n, 0.5, 2)

  L <- R <- numeric(n)
  for (i in seq_len(n)) {
    if      (X[i] <= Y[i]) { L[i] <- 0;    R[i] <- Y[i] }
    else if (X[i] >  Z[i]) { L[i] <- Z[i]; R[i] <- Inf  }
    else                    { L[i] <- Y[i]; R[i] <- Z[i] }
  }

  cat("\n===== H0 верна: X ~ Weibull(beta=2, eta=3) =====\n")
  res_h0 <- lb.gof.simple.weibull(L, R, beta = beta0, eta = eta0)

  cat("\n===== H1: неверный параметр формы, beta=1 =====\n")
  res_h1 <- lb.gof.simple.weibull(L, R, beta = 1, eta = eta0)

  # Weibull(1, eta) совпадает с Exp(eta) — проверка согласованности
  cat("\n===== Проверка: Exp(2) = Weibull(beta=1, eta=2) =====\n")
  X2 <- rexp(n, rate = 1/2)
  Y2 <- rexp(n, rate = 1/5)
  Z2 <- Y2 + rexp(n, rate = 1)
  L2 <- R2 <- numeric(n)
  for (i in seq_len(n)) {
    if      (X2[i] <= Y2[i]) { L2[i] <- 0;     R2[i] <- Y2[i] }
    else if (X2[i] >  Z2[i]) { L2[i] <- Z2[i]; R2[i] <- Inf   }
    else                      { L2[i] <- Y2[i]; R2[i] <- Z2[i] }
  }
  res_exp <- lb.gof.simple.weibull(L2, R2, beta = 1, eta = 2)
}
