# STGB33 · Labb 4 · Kapitel 9
# Building the Regression Model I: Model Selection and Validation
# Surgical Unit — fyra prediktorer, base R

# ------------------------------------------------------------
# 1. Läs tränings- och valideringsdata
# ------------------------------------------------------------
train_raw <- read.csv("surgical_unit_training.csv")
valid_raw <- read.csv("surgical_unit_validation.csv")

prepare_data <- function(d) {
  transform(d,
    x1 = blood_clotting,
    x2 = prognostic_index,
    x3 = enzyme_test,
    x4 = liver_test,
    y = survival_time,
    lny = ln_survival
  )
}

train <- prepare_data(train_raw)
valid <- prepare_data(valid_raw)

str(train)
head(train)

# ------------------------------------------------------------
# 2. Undersök modellformen med den ursprungliga responsen Y
#    Kutner §9.2, Figure 9.2a–b
# ------------------------------------------------------------
model_y <- lm(y ~ x1 + x2 + x3 + x4, data = train)
summary(model_y)

res_y <- residuals(model_y)
fit_y <- fitted(model_y)

par(mfrow = c(2, 2))
plot(fit_y, res_y, xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs fitted — Y")
abline(h = 0, lty = 2)
qqnorm(res_y, main = "Normal Q-Q plot — Y")
qqline(res_y)
plot(train$x3, res_y, xlab = "X3: enzyme function test",
     ylab = "Residuals", main = "Residuals vs X3 — Y")
abline(h = 0, lty = 2)
plot(seq_along(res_y), res_y, xlab = "Case", ylab = "Residuals",
     main = "Residuals by case — Y")
abline(h = 0, lty = 2)
par(mfrow = c(1, 1))

# Diagnostiken visar att Y-skalan inte är lämplig:
# spridningen ökar och residualerna avviker från normalitet.

# ------------------------------------------------------------
# 3. Transformera till ln(Y), skatta om och kontrollera igen
#    Kutner §9.2, Figure 9.2c–d och Figure 9.3
# ------------------------------------------------------------
train$lny <- log(train$y)
valid$lny <- log(valid$y)

round(cor(train[, c("lny", "x1", "x2", "x3", "x4")]), 4)

model_log_full <- lm(lny ~ x1 + x2 + x3 + x4, data = train)
summary(model_log_full)

res_log <- residuals(model_log_full)
fit_log <- fitted(model_log_full)

par(mfrow = c(2, 2))
plot(fit_log, res_log, xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs fitted — ln(Y)")
abline(h = 0, lty = 2)
qqnorm(res_log, main = "Normal Q-Q plot — ln(Y)")
qqline(res_log)
plot(train$x3, res_log, xlab = "X3: enzyme function test",
     ylab = "Residuals", main = "Residuals vs X3 — ln(Y)")
abline(h = 0, lty = 2)
plot(seq_along(res_log), res_log, xlab = "Case", ylab = "Residuals",
     main = "Residuals by case — ln(Y)")
abline(h = 0, lty = 2)
par(mfrow = c(1, 1))

# ------------------------------------------------------------
# 4. All possible subsets för X1–X4: 2^4 = 16 modeller
#    Kutner §9.3
# ------------------------------------------------------------
# Installera paketet en gång om det saknas:
# install.packages("olsrr")
library(olsrr)

# Endast en full modell behöver skattas.
model_full_subsets <- lm(lny ~ x1 + x2 + x3 + x4, data = train)

# Undersök alla möjliga kombinationer av prediktorerna.
subsets4 <- ols_step_all_possible(model_full_subsets)
subsets4

# Diagnosdiagram för modellvalskriterierna.
plot(subsets4)

# Paketet redovisar de 15 modeller som innehåller minst en prediktor.
# Interceptmodellen är den sextonde teoretiskt möjliga modellen.
model_intercept <- lm(lny ~ 1, data = train)
summary(model_intercept)

# Välj modellen med högst justerat R2 inom varje modellstorlek.
ordered_subsets4 <- subsets4[order(subsets4$n, -subsets4$adjr), ]
best_subsets4 <- ordered_subsets4[!duplicated(ordered_subsets4$n), ]
best_subsets4

# ------------------------------------------------------------
# 5. Forward och backward selection med AIC
#    Kutner §9.4
# ------------------------------------------------------------
null_model <- lm(lny ~ 1, data = train)
full_model <- lm(lny ~ x1 + x2 + x3 + x4, data = train)

forward_model <- step(
  null_model,
  scope = list(lower = formula(null_model), upper = formula(full_model)),
  direction = "forward",
  trace = TRUE
)

backward_model <- step(
  full_model,
  direction = "backward",
  trace = TRUE
)

formula(forward_model)
formula(backward_model)
summary(forward_model)
summary(backward_model)

# Båda metoderna väljer ln(Y) ~ X1 + X2 + X3.

# ------------------------------------------------------------
# 6. Validera de två bästa kandidaterna
#    Modell A: vald modell X1–X3
#    Modell B: full jämförelsemodell X1–X4
#    Kutner §9.6, MSPR ekvation (9.20)
# ------------------------------------------------------------
model_a_train <- lm(lny ~ x1 + x2 + x3, data = train)
model_b_train <- lm(lny ~ x1 + x2 + x3 + x4, data = train)

pred_a <- predict(model_a_train, newdata = valid)
pred_b <- predict(model_b_train, newdata = valid)

# MSPR beräknas direkt, utan egen function.
mspr_a <- mean((valid$lny - pred_a)^2)
mspr_b <- mean((valid$lny - pred_b)^2)

validation_result <- data.frame(
  model = c("Modell A: X1+X2+X3", "Modell B: X1+X2+X3+X4"),
  MSE_training = c(summary(model_a_train)$sigma^2,
                   summary(model_b_train)$sigma^2),
  MSPR_validation = c(mspr_a, mspr_b)
)
validation_result_round <- validation_result
validation_result_round[, 2:3] <- round(validation_result_round[, 2:3], 4)
validation_result_round

# R-skript för MSE/MSPR-grafen.
plot_values <- rbind(validation_result$MSE_training,
                     validation_result$MSPR_validation)
barplot(plot_values, beside = TRUE,
        names.arg = c("Modell A", "Modell B"),
        col = c("#13b8a6", "#092235"),
        ylab = "Mean squared error", ylim = c(0, 0.10))
legend("topright", legend = c("MSE, training", "MSPR, validation"),
       fill = c("#13b8a6", "#092235"), bty = "n")

# Skatta samma två modellformer separat i valideringsdata.
# Jämför hela summary(), inte bara MSE och MSPR.
model_a_valid <- lm(lny ~ x1 + x2 + x3, data = valid)
model_b_valid <- lm(lny ~ x1 + x2 + x3 + x4, data = valid)

summary(model_a_train)
summary(model_a_valid)
summary(model_b_train)
summary(model_b_valid)

# Slutsats: Modell A är enklare och har stabila koefficienter.
# Modell B får marginellt lägre MSPR, men X4 är inte signifikant i
# vare sig tränings- eller valideringsdata och förbättringen är liten.
