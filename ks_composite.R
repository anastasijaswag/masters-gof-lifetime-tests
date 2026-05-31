# ==============================================================================
# ks_composite.R
# Критерий Колмогорова–Смирнова для сложной гипотезы (полная выборка)
# Stephens (1977), D'Agostino & Stephens (1986)
# ==============================================================================

source("parameter_estimates_Exp.R")
source("parameter_estimates_Weibull.R")

# ------------------------------------------------------------------------------
# Основная функция: модифицированный критерий КС (сложная гипотеза)
#
#   H0: данные из экспоненциального или Вейбулла распределения
#       (параметры оцениваются по МПО из полной выборки)
#
# Аргументы:
#   x      — числовой вектор наблюдений (полная выборка, без цензуры)
#   distr  — "exponential" или "weibull"
#   alpha  — уровень значимости: 0.10, 0.05, 0.025 или 0.01
# ------------------------------------------------------------------------------

ks.test.composite <- function(x, distr = "exponential", alpha = 0.05) {

  if (!distr %in% c("exponential", "weibull")) {
    stop("distr должен быть 'exponential' или 'weibull'")
  }

  n  <- length(x)
  xs <- sort(x)

  # Оценки параметров МПО
  if (distr == "exponential") {
    eta  <- mle_exp_complete(xs)          # возвращает скаляр eta_hat
  } else {
    est  <- mle_weibull_complete(xs)      # возвращает список
    eta  <- est$eta_hat
    beta <- est$beta_hat
  }

  # Теоретическая функция распределения
  if (distr == "exponential") {
    F0 <- 1 - exp(-xs / eta)
  } else {
    F0 <- 1 - exp(-(xs / eta)^beta)
  }

  # Статистика КС
  i_vec   <- seq_len(n)
  D_plus  <- i_vec / n - F0
  D_minus <- F0 - (i_vec - 1) / n
  D       <- max(c(D_plus, D_minus))

  # Модификация Стивенса
  D_mod <- (D - 0.2 / n) * (sqrt(n) + 0.26 + 0.5 / sqrt(n))

  # Критические значения: D'Agostino & Stephens (1986)
  # Таблица 4.17 (экспонента), Таблица 4.2 (Вейбулл)
  if (distr == "exponential") {
    crit_table <- c("0.10" = 0.990, "0.05" = 1.094,
                    "0.025" = 1.190, "0.01" = 1.308)
  } else {
    crit_table <- c("0.10" = 0.803, "0.05" = 0.874,
                    "0.025" = 0.939, "0.01" = 1.007)
  }

  alpha_str <- as.character(alpha)
  if (!alpha_str %in% names(crit_table)) {
    stop("alpha должен быть одним из: 0.10, 0.05, 0.025, 0.01")
  }
  D_crit <- crit_table[alpha_str]

  reject <- (D_mod > D_crit)

  result <- list(
    test         = "Колмогорова–Смирнова",
    distribution = ifelse(distr == "exponential",
                          "экспоненциальное", "Вейбулла"),
    data_type    = "полная выборка",
    hypothesis   = ifelse(distr == "exponential",
                          "H0: выборка согласуется с экспоненциальным распределением",
                          "H0: выборка согласуется с распределением Вейбулла"),
    alpha        = alpha,
    n            = n,
    D            = D,
    D_mod        = D_mod,
    D_crit       = D_crit,
    reject       = reject
  )

  if (distr == "exponential") {
    result$eta_hat <- eta
  } else {
    result$eta_hat  <- eta
    result$beta_hat <- beta
  }

  cat("Критерий согласия:", result$test, "\n")
  cat("Распределение:", result$distribution, "\n")
  cat("Тип данных:", result$data_type, "\n")
  cat(result$hypothesis, "\n\n")

  cat("Оценки параметров:\n")
  cat("  eta_hat =", round(result$eta_hat, 4), "\n")
  if (distr == "weibull") {
    cat("  beta_hat =", round(result$beta_hat, 4), "\n")
  }
  cat("  alpha   =", round(result$alpha, 4), "\n\n")

  cat("Выборка:\n")
  cat("  n =", result$n, "\n\n")

  cat("Результат:\n")
  cat("  D*     =", round(result$D_mod, 4), "\n")
  cat("  D_crit =", round(result$D_crit, 4), "\n")
  cat("  Вывод ->",
      ifelse(result$reject, "H0 отклоняется", "H0 не отклоняется"), "\n\n")

  return(invisible(result))
}

# ==============================================================================
# Пример использования (не выполняется при source)
# ==============================================================================

if (FALSE) {

  set.seed(123)
  x_exp     <- rexp(200, rate = 1)
  x_weibull <- rweibull(200, shape = 2, scale = 3)

  cat("\n===== H0 верна: X ~ Exp(1), проверяем экспоненту =====\n")
  ks.test.composite(x_exp, distr = "exponential", alpha = 0.05)

  cat("\n===== H0 верна: X ~ Weibull(2, 3), проверяем Вейбулл =====\n")
  ks.test.composite(x_weibull, distr = "weibull", alpha = 0.05)

  cat("\n===== H1: X ~ Weibull(2, 3), проверяем экспоненту =====\n")
  ks.test.composite(x_weibull, distr = "exponential", alpha = 0.05)
}
