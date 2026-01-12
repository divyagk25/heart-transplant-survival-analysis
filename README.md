# Survival Analysis of Heart Transplant Patients
## An Applied Biostatistical Study Using Kaplan–Meier and Cox Models

---

## 1. Project Overview
This project investigates survival outcomes among heart transplant patients using classical survival analysis techniques. The primary objective is to assess whether receiving a heart transplant significantly improves patient survival after adjusting for relevant clinical covariates such as age, waiting time, and prior heart conditions.

The analysis applies Kaplan–Meier survival estimation, log-rank tests, and Cox proportional hazards models, along with rigorous diagnostic checks and model comparisons.

---

## 2. Research Questions
- Does receiving a heart transplant improve survival compared to not receiving one?
- How do age, waiting time, and prior medical conditions affect survival time?
- Are the effects of these covariates constant over time?
- Which survival model best explains the observed data?

---

## 3. Dataset Description
The analysis uses the `heart_transplant_tbl_df` dataset from the **MedDataSets** R package, containing data on **103 patients**.

Key variables include:
- `survtime`: survival time (days)
- `survived`: survival status (alive / dead)
- `transplant`: transplant status (control vs transplant)
- `age`: age at acceptance
- `wait`: waiting time before transplant
- `prior`: prior heart condition indicator

Approximately **73% of patients died** during the follow-up period, with substantial differences observed between transplant and non-transplant groups.

---

## 4. Data Preprocessing
- Missing values in waiting time were imputed using the median.
- Categorical variables were encoded numerically for survival modeling.
- Survival objects were constructed using time-to-event and censoring indicators.

These steps ensured consistency and compatibility with survival analysis methods.

---

## 5. Kaplan–Meier Survival Analysis
An overall Kaplan–Meier curve showed a steady decline in survival probability over time, with survival dropping sharply within the first 500 days.

Kaplan–Meier curves stratified by transplant group revealed:
- Significantly longer survival among transplant recipients
- A log-rank test p-value < 0.001, confirming statistically significant differences
- Median survival of approximately 21 days for controls versus 207 days for transplant patients

Log–log survival plots were approximately parallel, suggesting partial support for the proportional hazards assumption.

---

## 6. Cox Proportional Hazards Models

### 6.1 Simple Cox Model
A univariable Cox model including transplant status showed:
- Hazard Ratio ≈ 0.27
- Approximately 73% reduction in hazard of death for transplant patients

However, diagnostic tests indicated violation of the proportional hazards assumption.

---

### 6.2 Multivariable Cox Model
A multivariable model adjusting for age, waiting time, and prior condition showed:
- Transplant remained strongly protective
- Age was significantly associated with increased hazard
- PH assumption violations for transplant, wait time, and prior condition

---

### 6.3 Stratified Cox Model
To address non-proportional hazards, a stratified Cox model was fitted:
- Stratification on wait time and prior condition
- Transplant effect remained protective
- Age remained a significant predictor (HR ≈ 1.04 per year)

---

### 6.4 Interaction Model
An interaction between transplant status and age was tested:
- Interaction term was not statistically significant
- No evidence that transplant effect varies by age

---

## 7. Model Comparison
Models were compared using Akaike Information Criterion (AIC):

| Model | AIC |
|------|-----|
| Simple Cox | 572.3 |
| Multivariable Cox | 552.8 |
| Stratified Cox | **194.0** |
| Interaction Model | 195.5 |

The **stratified Cox model** provided the best fit.

---

## 8. Causal Framework
A Directed Acyclic Graph (DAG) was constructed to represent assumed causal relationships between transplant status, survival, and covariates. This helped clarify confounding structures and justified adjustment choices.

---

## 9. Key Conclusions
- Heart transplantation significantly improves survival
- Older age increases mortality risk
- PH assumption violations require stratified modeling
- No evidence of age–transplant interaction
- Stratified Cox model is most appropriate for this dataset

---

## 10. Limitations
- Small sample size (n = 103)
- Some covariates violate proportional hazards
- Observational data limits causal interpretation

---

## 11. Repository Structure
```text
code/        → R script with full survival analysis  
reports/     → Detailed statistical report  
images/      → Plots used in this README  
data/        → Data source description  
