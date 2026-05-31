# Критерии согласия для моделей времени жизни в условиях цензурирования в прикладных задачах

Репозиторий содержит реализацию на языке R критериев согласия для проверки экспоненциального распределения и распределения Вейбулла по данным с различными типами цензурирования.

**Выпускная квалификационная работа**  
Казакова А. — Санкт-Петербургский государственный университет, 2026

---

## Требования

- **R** версии 4.0 и выше
- Пакеты: [`interval`](https://cran.r-project.org/package=interval), [`goftest`](https://cran.r-project.org/package=goftest)

```r
install.packages(c("interval", "goftest"))
```

Файлы `right_censoring_Exp.R` и `right_censoring_Weibull.R` используют только базовые функции R и не требуют дополнительных пакетов.

---

## Структура репозитория

### Вспомогательные файлы

| Файл | Описание |
|------|----------|
| `utils_interval.R` | Общие утилиты для критериев при интервальной цензуре: работа с NPMLE, трапециевидное интегрирование, выбор N |
| `parameter_estimates_Exp.R` | МП-оценки параметра `eta` экспоненциального распределения: полная, правоцензурированная и интервально-цензурированная выборки |
| `parameter_estimates_Weibull.R` | МП-оценки параметров `(beta, eta)` распределения Вейбулла: полная, правоцензурированная и интервально-цензурированная выборки |

### Критерии для полной выборки

| Файл | Функция | Описание |
|------|---------|----------|
| `ks_composite.R` | `ks.test.composite()` | Критерий Колмогорова–Смирнова, сложная гипотеза (Stephens, 1977) |

### Критерии для интервально-цензурированных данных

| Файл | Функция | Гипотеза |
|------|---------|----------|
| `simple_test_Exp.R` | `lb.gof.simple.exp()` | Простая: F = Exp(eta), eta известно |
| `simple_test_Weibull.R` | `lb.gof.simple.weibull()` | Простая: F = Weibull(beta, eta), параметры известны |
| `complex_test_Exp.R` | `pf.lb.gof.exp()` | Сложная: F ∈ {Exp(eta) \| eta > 0} |
| `complex_test_Weibull.R` | `pf.lb.gof.weibull()` | Сложная: F ∈ {Weibull(beta, eta) \| beta, eta > 0} |

Все четыре критерия реализуют метод LB-GOF / PF-LB-GOF на основе непараметрической оценки максимального правдоподобия (NPMLE, алгоритм Тёрнбулла) и статистики Крамера–фон Мизеса (Ren, 2003).

### Критерии для правоцензурированных данных (II тип)

| Файл | Функция | Описание |
|------|---------|----------|
| `right_censoring_Exp.R` | `gof.exp.censored()` | Критерии W² и A² для экспоненциального распределения, табличные критические значения (D'Agostino & Stephens, 1986) |
| `right_censoring_Weibull.R` | `gof.weibull.censored()` | Параметрический бутстрэп-критерий для распределения Вейбулла (Zhu, 2020) |

---

## Быстрый старт

Перед запуском установите рабочую директорию в папку с файлами:

```r
setwd("путь/к/папке")
```

### Интервальная цензура — сложная гипотеза (Вейбулл)

```r
source("complex_test_Weibull.R")

# L, R — левые и правые границы интервалов наблюдения
# левая цензура: L = 0;  правая цензура: R = Inf
result <- pf.lb.gof.weibull(L, R, alpha = 0.05)
```

### Интервальная цензура — сложная гипотеза (экспонента)

```r
source("complex_test_Exp.R")

result <- pf.lb.gof.exp(L, R, alpha = 0.05)
```

### Интервальная цензура — простая гипотеза

```r
source("simple_test_Weibull.R")

# beta и eta — известные параметры нулевой гипотезы
result <- lb.gof.simple.weibull(L, R, beta = 2, eta = 3)
```

### Правая цензура II типа — Вейбулл (бутстрэп Жу)

```r
source("right_censoring_Weibull.R")

# x_obs — наблюдаемые отказы (первые m из n), n — полный объём
result <- gof.weibull.censored(x_obs, n = 100, B = 499)
```

### Правая цензура II типа — экспонента (Стивенс)

```r
source("right_censoring_Exp.R")

result <- gof.exp.censored(x_obs, n = 100, alpha = 0.05)
```

### Полная выборка — критерий КС

```r
source("ks_composite.R")   # требует parameter_estimates_Exp.R и parameter_estimates_Weibull.R

ks.test.composite(x, distr = "weibull", alpha = 0.05)
```

### МП-оценки параметров

```r
source("parameter_estimates_Exp.R")
source("parameter_estimates_Weibull.R")

# Полная выборка
eta <- mle_exp_complete(x)
mle_weibull_complete(x)       # возвращает список: eta_hat, beta_hat

# Правоцензурированная выборка
mle_exp_censored(times, indicators)
mle_weibull_censored(times, indicators)

# Интервально-цензурированная выборка
mle_exp_interval(t1, t2)
mle_weibull_interval(t1, t2)
```

---

## Параметризация распределения Вейбулла

Во всех файлах используется единая параметризация:

$$F(x;\, \beta, \eta) = 1 - \exp\\left(-\left(\frac{x}{\eta}\right)^{\\beta}\right), \quad x \geq 0$$

- `beta` — параметр формы
- `eta` — параметр масштаба

Совпадает с `pweibull(x, shape = beta, scale = eta)` в базовом R.

---

## Источники

- Ren, J.-J. (2003). Goodness of fit tests with interval censored data. *Scandinavian Journal of Statistics*, 30(1), 211–226.
- Zhu, T. (2020). Parametric bootstrap for goodness-of-fit testing with censored data. *Journal of Statistical Computation and Simulation*, 90(16), 2968–2983.
- D'Agostino, R. B., & Stephens, M. A. (1986). *Goodness-of-Fit Techniques*. Marcel Dekker.
- Turnbull, B. W. (1976). The empirical distribution function with arbitrarily grouped, censored and truncated data. *Journal of the Royal Statistical Society B*, 38(3), 290–295.
