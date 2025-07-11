# 1. Define the model parameters
ve_true = 0.90  # True vaccine efficacy for an all-or-nothing vaccine (90%)
force_of_infection = 0.05  # 5 cases per 100 person-years among susceptibles
population_size = 10000  # Size of each group
time_years = 1  # Duration of the study in years

# 2. Calculate person-time at risk for each group
person_time_u = population_size * time_years
# For a cohort study, all vaccinated individuals are considered in the denominator
# for person-time, even though a fraction is immune.
person_time_v = population_size * time_years

# 3. Calculate expected number of cases and Incidence Rates
# Unvaccinated Group: All are susceptible
cases_u = population_size * force_of_infection * time_years
ir_u = cases_u / person_time_u

# Vaccinated Group: Only the (1 - VE_true) fraction is susceptible
susceptible_v_fraction = 1 - ve_true
cases_v = (population_size * susceptible_v_fraction) * force_of_infection * time_years
ir_v = cases_v / person_time_v

# 4. Calculate IRR and the estimated VE from the IRR
irr = ir_v / ir_u
ve_estimated = 1 - irr

# 5. Print results and compare
print("--- Model Parameters ---")
print(f"True All-or-Nothing VE (VE_true): {ve_true:.2f}")
print(f"Force of Infection (among susceptibles): {force_of_infection}")
print("\n--- Calculated Incidence Rates ---")
print(f"Incidence Rate in Unvaccinated (IR_u): {ir_u:.4f}")
print(f"Incidence Rate in Vaccinated (IR_v):   {ir_v:.4f}")
print("\n--- VE Calculation ---")
print(f"Incidence Rate Ratio (IRR = IR_v / IR_u): {irr:.2f}")
# The user wants the full equation in the output.
print(f"Estimated VE = 1 - IRR")
print(f"             = 1 - {irr:.2f}")
print(f"             = {ve_estimated:.2f}")

print("\n--- Conclusion ---")
if abs(ve_estimated - ve_true) < 1e-9:
    print(f"The estimated VE ({ve_estimated:.2f}) correctly estimates the true VE ({ve_true:.2f}).")
else:
    print(f"The estimated VE ({ve_estimated:.2f}) does not correctly estimate the true VE ({ve_true:.2f}).")
