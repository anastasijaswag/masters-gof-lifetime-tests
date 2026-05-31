# ==============================================================================
# complex_test_Exp.R
# PF-LB-GOF критерий: сложная гипотеза H0: F in {Exp(eta) | eta > 0}
# Ren (2003), Scand. J. Statist., Section 4, теорема 3, ур. (19)–(22)
# ==============================================================================

source("utils_interval.R")

# ------------------------------------------------------------------------------
# Оценка параметра eta из NPMLE: theta_hat_n
# Remark 2 (с. 217): theta_hat_n — гладкий функционал от F_hat_n
# Для Exp(eta): eta = E[X], оценивается как среднее F_hat_n
# ------------------------------------------------------------------------------

estimate_theta_exp <- function(fit) {
  ap      <- get_npmle_atoms(fit)
  eta_hat <- sum(ap$atoms * ap$probs)   # среднее NPMLE
  eta_hat
}

# ------------------------------------------------------------------------------
# Статистика Крамера–фон Мизеса для сложной гипотезы (ур. (19))
#
#   T~*_m = m * int (F*_nm(x) - F0(x; eta_hat))^2 dF0(x; eta_hat)
#
# Через PIT: U_i = F0(X*_i; eta_hat),
#
#   T~*_m = 1/(12m) + sum [ U_{(i)} - (2i-1)/(2m) ]²
# ------------------------------------------------------------------------------

cvm_stat_composite_exp <- function(x_star, eta_hat) {
  m <- length(x_star)
  u <- sort(pexp(x_star, rate = 1 / eta_hat))
  1 / (12 * m) + sum((u - (2 * seq_len(m) - 1) / (2 * m))^2)
}

# ------------------------------------------------------------------------------
# Выбор m для сложной гипотезы (ур. (22))
#   m = max{ n^gamma, m_hat_tilde }
# где m_hat_tilde из (16)–(17) с заменой F0 на F0(·; theta_hat_n)
# ------------------------------------------------------------------------------

choose_m_composite_exp <- function(fit, n, eta_hat,
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

  x_max  <- qexp(1 - 1e-6, rate = 1 / eta_hat)
  x_grid <- seq(0, x_max, length.out = 5000)

  # F0 и f0 при оценённом theta_hat_n
  F0_vals <- pexp(x_grid, rate = 1 / eta_hat)
  f0_vals <- dexp(x_grid, rate = 1 / eta_hat)
  Fn_vals <- stepfun_npmle(fit, x_grid)

  # r_n = int (F_hat_n - F0(·; eta_hat))^2 dF0(·; eta_hat)  (ур. (17))
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
# Основная функция: PF-LB-GOF критерий (сложная гипотеза, экспонента)
#
#   H0: F in {Exp(eta) | eta > 0}  против  H1: F не принадлежит семейству
#
# Аргументы:
#   L, R  — левые и правые границы интервалов наблюдения;
#            левая цензура: L = 0;  правая цензура: R = Inf
# ------------------------------------------------------------------------------

pf_lb_gof_exp <- function(L, R,
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
  stopifnot(length(R) == n, all(L <= R))

  # Шаг 1: NPMLE
  fit <- icfit(L, R)

  # Шаг 2: оценка eta_hat из F_hat_n (Remark 2)
  eta_hat <- estimate_theta_exp(fit)

  # Шаг 3: выбор m с F0(·; eta_hat)  (ур. (22))
  m       <- choose_m_composite_exp(
    fit = fit, n = n, eta_hat = eta_hat,
    alpha = alpha, epsilon = epsilon,
    eta_param = eta_param, gamma = gamma,
    K_pilot = K_pilot
  )
  C_alpha <- qCvM(1 - alpha)

  # Шаг 4: оценка p~_n = P_n{ T~*_m >= C_alpha }  (ур. (11) + (19))
  T_pilot <- numeric(n_pilot)
  for (k in seq_len(n_pilot)) {
    x_star     <- sample_from_npmle(fit, m)
    T_pilot[k] <- cvm_stat_composite_exp(x_star, eta_hat)
  }
  p_n <- mean(T_pilot >= C_alpha)

  # Шаг 5: выбор N~ (ур. (12) с p~_n)
  N <- choose_N_lb(p_n = p_n, alpha = alpha, rho = rho,
                   N_max = N_max, N_min = N_min)$N

  # Шаг 6: внешний цикл PF-LB-GOF, вычисление W~_bar  (ур. (21))
  W_j <- integer(N)
  for (j in seq_len(N)) {
    x_star <- sample_from_npmle(fit, m)
    T_star <- cvm_stat_composite_exp(x_star, eta_hat)
    W_j[j] <- as.integer(T_star >= C_alpha)
  }
  W_bar <- mean(W_j)

  # Шаг 7: правило отклонения (ур. (21))
  z_alpha_rho <- qnorm(1 - (alpha - rho))
  threshold   <- alpha + z_alpha_rho * sqrt(alpha * (1 - alpha) / N)
  reject      <- (W_bar >= threshold)

  if (verbose) {
    cat("Критерий согласия: Крамера–фон Мизеса (PF-LB-GOF)\n")
    cat("Распределение: экспоненциальное\n")
    cat("Тип данных: интервально-цензурированная выборка\n")
    cat("H0: F ∈ {Exp(eta) | eta > 0}\n\n")

    cat("Оценки параметров:\n")
    cat("  eta_hat =", round(eta_hat, 4), "\n")
    cat("  alpha   =", round(alpha, 4), "\n\n")

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
    eta_hat   = eta_hat,
    npmle_fit = fit,
    mode      = "pf-lb-gof-exp"
  )
}

# ==============================================================================
# Пример использования (не выполняется при source)
# ==============================================================================

if (FALSE) {

  set.seed(2024)
  n <- 200

  make_intervals <- function(X, n) {
    Y <- rexp(n, rate = 1/5)
    Z <- Y + rexp(n, rate = 1)
    L <- R <- numeric(n)
    for (i in seq_len(n)) {
      if      (X[i] <= Y[i]) { L[i] <- 0;    R[i] <- Y[i] }
      else if (X[i] >  Z[i]) { L[i] <- Z[i]; R[i] <- Inf  }
      else                    { L[i] <- Y[i]; R[i] <- Z[i] }
    }
    list(L = L, R = R)
  }

  # H0 верна: X ~ Exp(2)
  obs <- make_intervals(rexp(n, rate = 1/2), n)
  cat("\n===== H0 верна: X ~ Exp(2) =====\n")
  pf_lb_gof_exp(obs$L, obs$R)

  # H1: X ~ Weibull(2, 3) — не экспонента
  obs2 <- make_intervals(rweibull(n, shape = 2, scale = 3), n)
  cat("\n===== H1: X ~ Weibull(2, 3) =====\n")
  pf_lb_gof_exp(obs2$L, obs2$R)
}
