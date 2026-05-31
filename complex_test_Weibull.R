# ==============================================================================
# complex_test_Weibull.R
# PF-LB-GOF критерий: сложная гипотеза H0: F in {Weibull(beta, eta) | beta, eta > 0}
# Ren (2003), Scand. J. Statist., Section 4, теорема 3, ур. (19)–(22)
#
# Параметризация:
#   F(x; beta, eta) = 1 - exp( -(x/eta)^beta ),   x >= 0
#   beta — параметр формы
#   eta  — параметр масштаба
# ==============================================================================

source("utils_interval.R")

# ------------------------------------------------------------------------------
# Оценка параметров Вейбулла из NPMLE: theta_hat_n = (beta_hat, eta_hat)
#
# Используется профильный МПО (Remark 2: theta_hat_n — гладкий функционал F_hat_n):
#   - для фиксированного beta параметр eta профилируется аналитически:
#       eta_hat(beta) = ( sum_j p_j * a_j^beta )^{1/beta}
#   - затем профильное правдоподобие максимизируется по beta (L-BFGS-B)
#   - при сбое оптимизатора — откат к методу моментов
#
# Метод моментов используется как начальное приближение и как резервный вариант:
#   CV^2 = Gamma(1+2/beta) / Gamma(1+1/beta)^2 - 1  =>  решаем для beta через uniroot
#
# Условие (AS2) из теоремы 3 выполнено: theta_hat_n — гладкий функционал от F_hat_n
# (Geskus & Groeneboom, 1999, теорема 3.2), плотность Вейбулла равномерно
# непрерывна по (beta, eta) на компактах.
# ------------------------------------------------------------------------------

estimate_theta_weibull_mom <- function(atoms, probs) {
  mu1     <- sum(atoms * probs)
  mu2     <- sum(atoms^2 * probs)
  var_hat <- mu2 - mu1^2

  if (var_hat <= 0 || mu1 <= 0) {
    return(list(beta = 1, eta = max(mu1, .Machine$double.eps)))
  }

  cv2 <- var_hat / mu1^2

  g <- function(b) gamma(1 + 2/b) / gamma(1 + 1/b)^2 - 1 - cv2

  root <- tryCatch(
    uniroot(g, interval = c(0.05, 200), tol = 1e-10),
    error = function(e) NULL
  )

  beta <- if (is.null(root)) 1 else root$root
  eta  <- mu1 / gamma(1 + 1/beta)

  list(beta = beta, eta = eta)
}

estimate_theta_weibull <- function(fit) {
  ap    <- get_npmle_atoms(fit)
  atoms <- ap$atoms
  probs <- ap$probs

  # Атомы в нуле дают log(0) = -Inf в правдоподобии — исключаем
  pos <- atoms > 0
  if (!all(pos)) {
    dropped_mass <- sum(probs[!pos])
    atoms <- atoms[pos]
    probs <- probs[pos]
    if (length(atoms) < 2 || dropped_mass > 0.5) {
      return(estimate_theta_weibull_mom(ap$atoms, ap$probs))
    }
    probs <- probs / sum(probs)
  }

  # Начальное приближение из метода моментов
  mom      <- estimate_theta_weibull_mom(atoms, probs)
  b_start  <- mom$beta

  # Профильное логарифмическое правдоподобие по beta
  # l_p(beta) = log beta - beta * log eta_hat(beta) + (beta-1) * sum p_j log a_j - 1
  sum_p_log_a <- sum(probs * log(atoms))

  neg_log_lik_profile <- function(b) {
    if (b <= 0 || !is.finite(b)) return(1e10)
    # sum p_j * a_j^beta — вычисляем устойчиво через log-sum-exp
    log_terms    <- log(probs) + b * log(atoms)
    max_lt       <- max(log_terms)
    log_sum_p_ab <- max_lt + log(sum(exp(log_terms - max_lt)))
    log_eta_hat  <- log_sum_p_ab / b
    ll <- log(b) - b * log_eta_hat + (b - 1) * sum_p_log_a - 1
    if (!is.finite(ll)) return(1e10)
    -ll
  }

  opt <- tryCatch(
    optim(
      par     = b_start,
      fn      = neg_log_lik_profile,
      method  = "L-BFGS-B",
      lower   = 0.05,
      upper   = 200,
      control = list(factr = 1e7)
    ),
    error = function(e) NULL
  )

  if (is.null(opt) || opt$convergence != 0 || !is.finite(opt$value)) {
    return(estimate_theta_weibull_mom(atoms, probs))
  }

  beta       <- opt$par
  log_terms  <- log(probs) + beta * log(atoms)
  max_lt     <- max(log_terms)
  log_sum_p_ab <- max_lt + log(sum(exp(log_terms - max_lt)))
  eta        <- exp(log_sum_p_ab / beta)

  if (!is.finite(beta) || !is.finite(eta) || beta <= 0 || eta <= 0) {
    return(estimate_theta_weibull_mom(atoms, probs))
  }

  list(beta = beta, eta = eta)
}

# ------------------------------------------------------------------------------
# Статистика Крамера–фон Мизеса для сложной гипотезы (ур. (19))
#
#   T~*_m = m * int (F*_nm(x) - F0(x; theta_hat))^2 dF0(x; theta_hat)
#
# Через PIT: U_i = F0(X*_i; theta_hat),
#
#   T~*_m = 1/(12m) + sum [ U_{(i)} - (2i-1)/(2m) ]²
# ------------------------------------------------------------------------------

cvm_stat_composite_weibull <- function(x_star, beta_hat, eta_hat) {
  m <- length(x_star)
  u <- sort(pweibull(x_star, shape = beta_hat, scale = eta_hat))
  1 / (12 * m) + sum((u - (2 * seq_len(m) - 1) / (2 * m))^2)
}

# ------------------------------------------------------------------------------
# Выбор m для сложной гипотезы (ур. (22))
#   m = max{ n^gamma, m_hat_tilde }
# где m_hat_tilde из (16)–(17) с заменой F0 на F0(·; theta_hat_n)
# ------------------------------------------------------------------------------

choose_m_composite_weibull <- function(fit, n, beta_hat, eta_hat,
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

  x_max  <- qweibull(1 - 1e-6, shape = beta_hat, scale = eta_hat)
  # Сетка начинается от machine epsilon, чтобы избежать f(0) = Inf при beta < 1
  x_grid <- seq(.Machine$double.eps, x_max, length.out = 5000)

  # F0 и f0 при оценённом theta_hat_n
  F0_vals <- pweibull(x_grid, shape = beta_hat, scale = eta_hat)
  f0_vals <- dweibull(x_grid, shape = beta_hat, scale = eta_hat)
  Fn_vals <- stepfun_npmle(fit, x_grid)

  # r_n = int (F_hat_n - F0(·; theta_hat))^2 dF0(·; theta_hat)  (ур. (17))
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
# Основная функция: PF-LB-GOF критерий (сложная гипотеза, Вейбулл)
#
#   H0: F in {Weibull(beta, eta) | beta, eta > 0}  против  H1: F не принадлежит семейству
#
# Аргументы:
#   L, R  — левые и правые границы интервалов наблюдения;
#            левая цензура: L = 0;  правая цензура: R = Inf
# ------------------------------------------------------------------------------

pf.lb.gof.weibull <- function(L, R,
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

  # Шаг 2: оценка theta_hat_n = (beta_hat, eta_hat) из F_hat_n (Remark 2)
  theta_hat <- estimate_theta_weibull(fit)
  beta_hat  <- theta_hat$beta
  eta_hat   <- theta_hat$eta

  # Шаг 3: выбор m с F0(·; theta_hat_n)  (ур. (22))
  m       <- choose_m_composite_weibull(
    fit = fit, n = n,
    beta_hat = beta_hat, eta_hat = eta_hat,
    alpha = alpha, epsilon = epsilon,
    eta_param = eta_param, gamma = gamma,
    K_pilot = K_pilot
  )
  C_alpha <- qCvM(1 - alpha)

  # Шаг 4: оценка p~_n = P_n{ T~*_m >= C_alpha }  (ур. (11) + (19))
  T_pilot <- numeric(n_pilot)
  for (k in seq_len(n_pilot)) {
    x_star     <- sample_from_npmle(fit, m)
    T_pilot[k] <- cvm_stat_composite_weibull(x_star, beta_hat, eta_hat)
  }
  p_n <- mean(T_pilot >= C_alpha)

  # Шаг 5: выбор N~ (ур. (12) с p~_n)
  N <- choose_N_lb(p_n = p_n, alpha = alpha, rho = rho,
                   N_max = N_max, N_min = N_min)$N

  # Шаг 6: внешний цикл PF-LB-GOF, вычисление W~_bar  (ур. (21))
  W_j <- integer(N)
  for (j in seq_len(N)) {
    x_star  <- sample_from_npmle(fit, m)
    T_star  <- cvm_stat_composite_weibull(x_star, beta_hat, eta_hat)
    W_j[j]  <- as.integer(T_star >= C_alpha)
  }
  W_bar <- mean(W_j)

  # Шаг 7: правило отклонения (ур. (21))
  z_alpha_rho <- qnorm(1 - (alpha - rho))
  threshold   <- alpha + z_alpha_rho * sqrt(alpha * (1 - alpha) / N)
  reject      <- (W_bar >= threshold)

  if (verbose) {
    cat("Критерий согласия: Крамера–фон Мизеса (PF-LB-GOF)\n")
    cat("Распределение: Вейбулла\n")
    cat("Тип данных: интервально-цензурированная выборка\n")
    cat("H0: F ∈ {Weibull(beta, eta) | beta, eta > 0}\n\n")

    cat("Оценки параметров:\n")
    cat("  beta_hat =", round(beta_hat, 4), "\n")
    cat("  eta_hat  =", round(eta_hat, 4), "\n")
    cat("  alpha    =", round(alpha, 4), "\n\n")

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
    beta_hat  = beta_hat,
    eta_hat   = eta_hat,
    npmle_fit = fit,
    mode      = "pf-lb-gof-weibull"
  )
}

# ==============================================================================
# Пример использования (не выполняется при source)
# ==============================================================================

if (FALSE) {

  set.seed(2024)
  n <- 200

  make_intervals <- function(X, n) {
    Y <- runif(n, 0, 6)
    Z <- Y + runif(n, 0.5, 2)
    L <- R <- numeric(n)
    for (i in seq_len(n)) {
      if      (X[i] <= Y[i]) { L[i] <- 0;    R[i] <- Y[i] }
      else if (X[i] >  Z[i]) { L[i] <- Z[i]; R[i] <- Inf  }
      else                    { L[i] <- Y[i]; R[i] <- Z[i] }
    }
    list(L = L, R = R)
  }

  # H0 верна: X ~ Weibull(beta=2, eta=3)
  obs <- make_intervals(rweibull(n, shape = 2, scale = 3), n)
  cat("\n===== H0 верна: X ~ Weibull(beta=2, eta=3) =====\n")
  pf.lb.gof.weibull(obs$L, obs$R)

  # H1: X ~ LogNormal(1, 0.5)
  obs2 <- make_intervals(rlnorm(n, meanlog = 1, sdlog = 0.5), n)
  cat("\n===== H1: X ~ LogNormal(1, 0.5) =====\n")
  pf.lb.gof.weibull(obs2$L, obs2$R)

  # Проверка согласованности: Exp(2) = Weibull(beta=1, eta=2)
  obs3 <- make_intervals(rexp(n, rate = 1/2), n)
  cat("\n===== Проверка: X ~ Exp(2) = Weibull(beta=1, eta=2) =====\n")
  res <- pf.lb.gof.weibull(obs3$L, obs3$R)
  cat("  beta_hat (ожидается ~1):", round(res$beta_hat, 4), "\n")
  cat("  eta_hat  (ожидается ~2):", round(res$eta_hat, 4), "\n")
}
