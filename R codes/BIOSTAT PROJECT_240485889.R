# Install and load necessary packages (run once)
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")
if (!require("MedDataSets")) install.packages("MedDataSets")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("gtsummary")) install.packages("gtsummary")
if (!require("ggsurvfit")) install.packages("ggsurvfit")

library(survival)
library(survminer)
library(MedDataSets)
library(dplyr)
library(ggplot2)
library(gtsummary)
library(ggsurvfit)

# Load the heart transplant dataset
data(heart_transplant_tbl_df)
heartdf <- as.data.frame(heart_transplant_tbl_df)

# View first few rows
head(heartdf)

# Check for missing values
colSums(is.na(heartdf))

# Impute missing 'wait' times with median 
median_wait <- median(heartdf$wait, na.rm = TRUE)
heartdf$wait[is.na(heartdf$wait)] <- median_wait

# Encode variables for survival analysis: convert categorical variables to numeric 
# survived: "alive" -> 1, "dead" -> 0 (convert to numeric)
heartdf$survived <- ifelse(heartdf$survived == "alive", 1, 0)
# prior: "yes" -> 1, "no" -> 0
heartdf$prior <- ifelse(heartdf$prior == "yes", 1, 0)
# Convert survival time to numeric
heartdf$survtime <- as.numeric(as.character(heartdf$survtime))


# Overall survival statistics
overall_stats <- heartdf %>%
  summarise(
    total_patients = n(),
    total_deaths = sum(survived == 0),
    death_rate = mean(survived == 0),
    median_surv_all = median(survtime)
  )
print(overall_stats)


# Create survival object which combines survival time and event status: event = death (1 - survived)
surv_obj <- Surv(time = heartdf$survtime, event = 1 - heartdf$survived)


# 1. Overall Kaplan-Meier survival curve: estimate survival probabilities over time
km_fit_overall <- survfit(surv_obj ~ 1)
ggsurvfit(km_fit_overall) +
  add_confidence_interval() +
  add_risktable() +
  labs(x = "Time (Days)",
       y = "Overall Survival Probability",
       title = "Overall Kaplan-Meier Survival Estimate")
#The overall Kaplan-Meier curve shows that patient survival declines steadily over time.



# Statistics by transplant group
group_stats <- heartdf %>%
  group_by(transplant) %>%
  summarise(
    n = n(),
    deaths = sum(survived == 0),
    death_rate = mean(survived == 0),
    median_surv = median(survtime),
    mean_age = mean(age),
    mean_wait = mean(wait)
  )
print(group_stats)


# 2. Kaplan-Meier survival curve stratified by transplant group: compare survival between patients who got a transplant vs those who didn't
km_fit_transplant <- survfit(surv_obj ~ transplant, data = heartdf)
ggsurvplot(km_fit_transplant,
           data = heartdf,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           title = "Kaplan-Meier Survival Curve by Transplant Group",
           xlab = "Time (Days)",
           ylab = "Survival Probability",
           palette = "Dark2",
           legend.title = "Transplant Group")

# Log-rank test to compare survival between transplant groups
log_rank_test <- survdiff(surv_obj ~ transplant, data = heartdf)
print(log_rank_test)
cat("Log-rank test p-value:", 1 - pchisq(log_rank_test$chisq, df = length(log_rank_test$n) - 1), "\n\n")

# Interpretation:
# The Kaplan-Meier plot shows that patients who received transplants generally survive longer.
# The log-rank test p-value < 0.001 indicates significant difference in survival and confirms that the survival difference is statistically significant.

#Log-Log Survival Plot
ggsurvplot(
  survfit(surv_obj ~ transplant, data = heartdf),
  fun = "cloglog",
  palette = c("blue", "red"),
  legend.title = "Transplant Type",
  ggtheme = theme_minimal(),
  title = "Log-Log Survival Plot by Transplant Group"
)
#Overall, this plot is approximately parallel and supports using a Cox proportional hazards model to analyze the survival difference between transplant groups.

# 3. Cox Proportional Hazards Model:  transplant only
cox_simple <- coxph(surv_obj ~ transplant, data = heartdf)
summary(cox_simple)  #model summary with hazard ratios and significance levels for transplant vs no transplant group
cat("Hazard ratio (HR) for transplant vs control:", exp(coef(cox_simple)), "\n\n")
# Test proportional hazards (PH) assumption of the Cox model
cox_test <- cox.zph(cox_simple)
print(cox_test)
# Plot Schoenfeld residuals for visual assessment
ggcoxzph(cox_test)

# Interpretation:
# Patients who received transplant had significantly lower hazard of death (HR << 1).
#The low p-value (0.002) suggests violation of the proportional hazards assumption for the transplant variable. 


# 4. Multivariable Cox Model: adjusting for age, wait time, prior condition to study how multiple variables together affect survival. 
cox_multi <- coxph(surv_obj ~ transplant + age + wait + prior, data = heartdf)
summary(cox_multi) #model summary with hazard ratios and significance levels

# Test proportional hazards (PH) assumption of the Cox model
ph_test <- cox.zph(cox_multi)
print(ph_test)
ggcoxzph(ph_test)

#  PH assumption violated (p < 0.05) for some variables. Hence, stratify by violating variables
# Here, transplant, wait, prior show violation

#stratifying:
cox_strat <- coxph(surv_obj ~ transplant + age + strata(prior) + strata(wait), data = heartdf)
summary(cox_strat)
ph_test_strat <- cox.zph(cox_strat) #check PH assumption
ggcoxzph(ph_test_strat)

# Interpretation:
# Stratified Cox model adjusts for non-proportional hazards in prior and wait.
# Having a heart transplant lowers the risk of death, but age increases hazard.

# 5. Add interaction between transplant and age:to see if transplant effect varies by age
cox_interaction <- coxph(surv_obj ~ transplant * age + strata(prior) + strata(wait), data = heartdf)
summary(cox_interaction)
ph_interaction <- cox.zph(cox_interaction)
ggcoxzph(ph_interaction)
#The interaction term is not statistically significant (p = 0.49), suggesting no evidence that the transplant effect differs by age. Age was significantly associated with survival (HR = 1.049, p = 0.015), indicating higher risk of death with increasing age.

# Fit Kaplan-Meier curves by prior condition
km_prior <- survfit(surv_obj ~ prior, data = heartdf)

# Plot survival curves by prior condition
ggsurvplot(km_prior, data = heartdf,
           title = "Survival by Prior Condition",
           xlab = "Time (Days)",
           ylab = "Survival Probability",
           risk.table = TRUE,
           pval = TRUE)
#p value=0.01
#Prior condition affects survival time: patients with prior conditions tend to have different survival outcomes compared to those without. 


# Compare models by AIC
AIC(cox_simple, cox_multi, cox_strat, cox_interaction)

#DAG
install.packages("dagitty")
install.packages("ggdag")
library(dagitty)
library(ggdag)

dag_simple <- dagify(
  survival ~ wait + transplant + age + prior,
  wait ~ transplant + age + prior,
  transplant ~ age + prior,
  exposure = "transplant",
  outcome = "survival"
)

ggdag(dag_simple, layout = "circle") +
  theme_minimal() +
  ggtitle("Simplified DAG: Heart Transplant → Survival")

