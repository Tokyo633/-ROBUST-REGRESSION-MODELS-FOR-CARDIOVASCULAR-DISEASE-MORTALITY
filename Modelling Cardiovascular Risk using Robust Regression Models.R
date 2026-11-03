getwd()
setwd("C:/Users/grace/OneDrive/Documents/School Stuff/Semester 4/Advanced Regression/Lab 4 assignment")
cvd<- read.csv("cvd_mortality.csv")
dim(cvd)
# ---- Setup ----
library(ggplot2)
library(GGally)    
library(car)        

cvd$community <- NULL        
cvd$id <- seq_len(nrow(cvd))

str(cvd)
summary(cvd)
colSums(is.na(cvd))

# ---- 1. Summary table (mean, sd, median, MAD, IQR) ----
vars <- c("rate", "pm25", "smoking", "income", "age65")
eda_tab <- t(sapply(cvd[vars], function(x)
  c(n = length(x), mean = mean(x), sd = sd(x), median = median(x),
    MAD = mad(x), IQR = IQR(x), min = min(x), max = max(x),
    skew = mean((x - mean(x))^3) / sd(x)^3)))
round(eda_tab, 2)

# ---- 2. Distributions: histogram + boxplot for each variable ----
op <- par(mfrow = c(2, 5), mar = c(4, 4, 2, 1))
for (v in vars) hist(cvd[[v]], main = v, xlab = v, col = "grey80", breaks = 20)
for (v in vars) boxplot(cvd[[v]], main = v, col = "grey80", horizontal = FALSE)
par(op)

# Which values are flagged by the boxplot rule (1.5 * IQR)?
for (v in vars) {
  out <- boxplot.stats(cvd[[v]])$out
  cat(v, ": ", length(out), "boxplot outliers ->", sort(out), "\n")
}

# Is income top-coded?
sum(cvd$income == max(cvd$income))
sort(cvd$income, decreasing = TRUE)[1:10]

# ---- 3. Robust z-scores: rows that are extreme on any variable ----
rz <- sapply(cvd[vars], function(x) (x - median(x)) / mad(x))
flag <- which(apply(abs(rz), 1, max) > 3)
cbind(id = flag, round(cvd[flag, vars], 1), round(rz[flag, ], 1))

# ---- 4. Relationships: scatterplot matrix with correlations ----
ggpairs(cvd[vars],
        lower = list(continuous = wrap("smooth", alpha = 0.4, size = 1)))
round(cor(cvd[vars]), 2)
round(cor(cvd[vars], method = "spearman"), 2)   # compare: big gaps hint at outlier influence

# ---- 5. Main relationship: rate vs pm25 (OLS line vs lowess vs robust-ish) ----
ggplot(cvd, aes(pm25, rate)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, colour = "red") +
  geom_smooth(method = "loess", se = FALSE, colour = "blue") +
  labs(title = "CVD mortality vs PM2.5",
       subtitle = "red = OLS, blue = loess",
       x = "PM2.5 (µg/m³)", y = "CVD deaths per 100,000")

# Same plot for each covariate in one panel
long <- reshape(cvd[c("id", "rate", "pm25", "smoking", "income", "age65")],
                varying = c("pm25", "smoking", "income", "age65"),
                v.names = "x", timevar = "var",
                times = c("pm25", "smoking", "income", "age65"),
                direction = "long")
ggplot(long, aes(x, rate)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", se = FALSE, colour = "red") +
  geom_smooth(method = "loess", se = FALSE, colour = "blue") +
  facet_wrap(~ var, scales = "free_x")

# ---- 6. Skew check: does a log transform help? (pm25, income, rate) ----
op <- par(mfrow = c(2, 3))
for (v in c("rate", "pm25", "income")) {
  hist(cvd[[v]], main = v, xlab = v, col = "grey80")
  hist(log(cvd[[v]]), main = paste("log", v), xlab = paste("log", v), col = "grey80")
}
par(op)

# ---- 7. Multicollinearity between the predictors ----
vif(lm(rate ~ pm25 + smoking + income + age65, data = cvd))

# ---- 8. Quick look at the suspicious pm25 points ----
cvd[order(-cvd$pm25)[1:8], ]
cvd[order(-cvd$rate)[1:8], ]

library(MASS)       
library(car)        
library(lmtest)     
library(sandwich)   
library(ggplot2)

cvd <- read.csv("cvd_mortality.csv")
cvd$id <- seq_len(nrow(cvd))
f <- rate ~ pm25 + smoking + income + age65

# ---- 1. Fit OLS ----
ols <- lm(f, data = cvd)
summary(ols)
confint(ols)

# Unadjusted vs adjusted PM2.5 slope (how much confounding there is)
cat("unadjusted pm25 slope:", coef(lm(rate ~ pm25, data = cvd))["pm25"], "\n")
cat("adjusted  pm25 slope:", coef(ols)["pm25"], "\n")

# ---- 2. Standard four-panel diagnostics ----
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1)); plot(ols, las = 1); par(mfrow = c(1, 1))

# ---- 3. Influence measures in one table ----
library(MASS)
n <- nrow(cvd); p <- length(coef(ols))
diag_tab <- data.frame(
  id       = cvd$id,
  rate     = cvd$rate, pm25 = cvd$pm25, smoking = cvd$smoking,
  income   = cvd$income, age65 = cvd$age65,
  leverage = hatvalues(ols),
  stud_res = studres(ols),
  cooks    = cooks.distance(ols),
  dffits   = dffits(ols)
)

lev_cut <- 2 * p / n          # leverage cutoff
cd_cut  <- 4 / n              # Cook's D cutoff
cat("Leverage cutoff:", round(lev_cut, 3),
    " Cook's D cutoff:", round(cd_cut, 4), "\n")

flag <- subset(diag_tab, leverage > lev_cut | abs(stud_res) > 3 | cooks > cd_cut)
flag[order(-flag$cooks), ]

# ---- 4. Residuals vs leverage, labelled + Cook's D index plot ----
ggplot(diag_tab, aes(leverage, stud_res)) +
  geom_point(aes(size = cooks), alpha = 0.5) +
  geom_hline(yintercept = c(-3, 3), linetype = 2) +
  geom_vline(xintercept = lev_cut, linetype = 2) +
  geom_text(data = flag, aes(label = id), vjust = -1, size = 3) +
  labs(x = "Leverage", y = "Studentized residual", size = "Cook's D")

plot(diag_tab$cooks, type = "h", ylab = "Cook's D", xlab = "Community id")
abline(h = cd_cut, col = "red", lty = 2)
text(flag$id, flag$cooks, labels = flag$id, pos = 3, cex = 0.7)



# ---- 6. Assumption tests ----
bptest(ols)                    # heteroscedasticity
shapiro.test(resid(ols))       # normality of residuals
vif(ols)                       # collinearity
car::ncvTest(ols)



# ---- 8. Sandwich (HC3) standard errors for the OLS fit ----
coeftest(ols, vcov = vcovHC(ols, type = "HC3"))

# ---- 9. Sensitivity: refit without the five high-PM2.5 communities ----
hi_pm <- cvd$pm25 > 20
ols_trim <- lm(f, data = cvd[!hi_pm, ])
summary(ols_trim)

## Huber model
huber_model <- rlm(
  rate ~ pm25 + smoking + income + age65,
  data = cvd,
  psi = psi.huber
)

summary(huber_model)

##Quantile regression
library(quantreg)

rq_model <- rq(
  rate ~ pm25 + smoking + income + age65,
  tau = 0.5,
  data = cvd
)

summary(rq_model)

## Thiel-Sen Estimation
install.packages("mblm")
library(mblm)

theilsen <- mblm(rate ~ pm25, data = cvd)

summary(theilsen)

bisquare_model <- rlm(
  rate ~ pm25 + smoking + income + age65,
  data = cvd,
  psi = psi.bisquare
)

summary(bisquare_model)
