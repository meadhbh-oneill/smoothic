# Dataset -----------------------------------------------------------------
library(tidyverse)



library(mlbench)
data("BostonHousing2")

raw_data <- BostonHousing2

dataset_raw <- raw_data %>%
  select(
    -town,
    -tract,
    -lon,
    -lat,
    -medv,
    -b
  ) %>%
  mutate(
    lcmedv = log(cmedv),
    lnox = log(nox),
    ldis = log(dis),
    ltax = log(tax),
    llstat = log(lstat),
    chast = as.numeric(as.character(chas)), # convert to numeric
    # chast = as.numeric(as.character(recode_factor(chas, "0" = "-1", "1" = "1"))),
    .keep = "unused"
  ) %>%
  relocate(lcmedv,
    .after = everything()
  )

# * Raw data ----
y_raw <- dataset_raw %>% pull(lcmedv)
x_raw <- dataset_raw %>%
  select(-lcmedv) %>%
  as.matrix()
x_raw_int <- as_tibble(x_raw) %>%
  add_column(
    inter = 1,
    .before = 1
  ) %>%
  as.matrix()

data <- cbind(x_raw, y_raw) %>% as.data.frame()

fit1 <- smoothic(formula = y_raw ~ .,
                data = data,
                model = "mpr",
                lambda = "log(n)",
                epsilon_1 = 0.3,
                epsilon_T = 1e-4,
                steps_T = 100,
                zero_tol = 1e-6,
                max_it = 1e4,
                family = "sgnd",
                optimizer = "nlm")

out


fit_sniffer <- smoothic(formula = y ~ .,
                        data = sniffer,
                        model = "mpr",
                        lambda = "log(n)",
                        epsilon_1 = 10,
                        epsilon_T = 1e-5,
                        steps_T = 100,
                        zero_tol = 1e-8,
                        max_it = 1e4,
                        family = "sgnd",
                        optimizer = "manual",
                        tau = 0.001,
                        tol = 1e-8,
                        initial_step = 10,
                        max_step_it = 1000,
                        kappa = 1.5)

fit_hprice <- smoothic(formula = lprice ~ .,
                        data = bostonhouseprice,
                        model = "mpr",
                        lambda = "log(n)",
                        epsilon_1 = 10,
                        epsilon_T = 1e-5,
                        steps_T = 100,
                        zero_tol = 1e-8,
                        max_it = 1e4,
                        family = "normal",
                        optimizer = "manual",
                        tau = 0.01,
                        tol = 1e-8,
                        initial_step = 10,
                        max_step_it = 1000)
