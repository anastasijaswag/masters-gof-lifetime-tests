# ==============================================================================
# utils_interval.R
# Общие вспомогательные функции для LB-GOF критериев (интервальная цензура)
# Используется: simple_test_Exp.R, simple_test_Weibull.R,
#               complex_test_Exp.R, complex_test_Weibull.R
# ==============================================================================

library(interval)   # icfit — NPMLE при интервальной цензуре
library(goftest)    # qCvM  — квантили распределения Крамера–фон Мизеса

# ------------------------------------------------------------------------------
# Работа с NPMLE (объект icfit)
# ------------------------------------------------------------------------------

# Извлекает опорные точки (атомы) и нормированные вероятностные массы из icfit
get_npmle_atoms <- function(fit) {
  probs     <- fit$pf
  left_pts  <- fit$intmap[1, ]
  right_pts <- fit$intmap[2, ]

  idx       <- is.finite(right_pts) & probs > 0
  probs     <- probs[idx]
  left_pts  <- left_pts[idx]
  right_pts <- right_pts[idx]

  if (length(probs) == 0L) {
    stop("NPMLE не имеет положительных масс на конечных точках.")
  }

  # Масса ставится в правый конец интервала Тёрнбулла
  atoms <- ifelse(left_pts == right_pts, left_pts, right_pts)

  ord   <- order(atoms)
  atoms <- atoms[ord]
  probs <- probs[ord]

  probs_by_atom <- tapply(probs, atoms, sum)
  atoms <- as.numeric(names(probs_by_atom))
  probs <- as.numeric(probs_by_atom)
  probs <- probs / sum(probs)

  list(atoms = atoms, probs = probs)
}

# Генерирует выборку размера m из дискретного NPMLE
sample_from_npmle <- function(fit, m) {
  ap <- get_npmle_atoms(fit)
  sample(ap$atoms, size = m, replace = TRUE, prob = ap$probs)
}

# Возвращает значения функции распределения NPMLE в точках x_grid
stepfun_npmle <- function(fit, x_grid) {
  ap <- get_npmle_atoms(fit)
  sf <- stepfun(ap$atoms, c(0, cumsum(ap$probs)), right = TRUE)
  sf(x_grid)
}

# ------------------------------------------------------------------------------
# Численное интегрирование методом трапеций
# ------------------------------------------------------------------------------

trapez <- function(x, y) {
  n <- length(x)
  if (n < 2L) return(0)
  dx <- diff(x)
  sum((y[-n] + y[-1]) * dx / 2)
}

# ------------------------------------------------------------------------------
# Выбор числа повторений N (Ren, 2003, ур. (12))
# Одинаков для простой и сложной гипотезы
# ------------------------------------------------------------------------------

choose_N_lb <- function(p_n,
                        alpha = 0.05,
                        rho   = 0.025,
                        N_max = 100000,
                        N_min = 30) {

  z_alpha     <- qnorm(1 - alpha)
  z_alpha_rho <- qnorm(1 - (alpha - rho))
  delta_z     <- z_alpha_rho - z_alpha

  denom <- (alpha - p_n)^2

  if (denom < 1e-10) {
    N <- N_max
  } else {
    N <- ceiling(p_n * (1 - p_n) * (delta_z^2) / denom)
    N <- max(N_min, min(N, N_max))
  }

  list(N = N)
}
