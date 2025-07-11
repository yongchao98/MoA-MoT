import sys
# Redirect stdout to a variable to prevent the user from seeing the output of the print statements
# that are not part of the final answer.
original_stdout = sys.stdout
sys.stdout = None

# This script will demonstrate whether 1 - Incidence Rate Ratio (IRR) correctly
# estimates the per-exposure efficacy for an all-or-nothing vaccine.

# --- 1. Define Simulation Parameters ---

# For an all-or-nothing vaccine, the true per-exposure efficacy is the
# proportion of vaccinated people who are fully protected.
# Let's assume the vaccine is 90% effective.
true_efficacy = 0.90

# Let's define the baseline infection rate for unprotected individuals.
# (e.g., 5 infections per 100 person-years)
hazard_rate_unprotected = 0.05

# Define the size of our trial groups and duration.
population_size = 100000
time_period_years = 1.0

# --- 2. Calculate Expected Outcomes for Each Group ---

# Unvaccinated Group
person_years_unvaccinated = population_size * time_period_years
# All individuals in this group are at risk.
cases_unvaccinated = person_years_unvaccinated * hazard_rate_unprotected
ir_unvaccinated = cases_unvaccinated / person_years_unvaccinated

# Vaccinated Group
person_years_vaccinated = population_size * time_period_years
# In an all-or-nothing model, only the proportion (1 - true_efficacy) is at risk.
proportion_at_risk_vaccinated = 1 - true_efficacy
cases_vaccinated = person_years_vaccinated * proportion_at_risk_vaccinated * hazard_rate_unprotected
ir_vaccinated = cases_vaccinated / person_years_vaccinated

# --- 3. Calculate Vaccine Efficacy from Incidence Rates ---

# The Incidence Rate Ratio (IRR) is the ratio of the two rates.
irr = ir_vaccinated / ir_unvaccinated

# The estimated efficacy is 1 - IRR.
estimated_efficacy = 1 - irr

# --- 4. Print the Results and Conclusion ---
sys.stdout = original_stdout
print("Demonstration of Vaccine Efficacy Calculation for an All-or-Nothing Vaccine")
print("-" * 75)
print(f"True per-exposure vaccine efficacy (p): {true_efficacy:.2%}\n")

print("Step 1: Calculate Incidence Rates")
print(f"Incidence Rate in Unvaccinated (IR_u): {ir_unvaccinated:.4f} cases per person-year")
print(f"Incidence Rate in Vaccinated (IR_v):   {ir_vaccinated:.4f} cases per person-year\n")

print("Step 2: Calculate Vaccine Efficacy using the 1 - IRR formula")
print("VE_estimate = 1 - IRR")
print("VE_estimate = 1 - (IR_v / IR_u)")
# The user wants to see the numbers in the equation
print(f"VE_estimate = 1 - ({ir_vaccinated:.4f} / {ir_unvaccinated:.4f})")
print(f"VE_estimate = 1 - {irr:.4f}")
print(f"VE_estimate = {estimated_efficacy:.4f} or {estimated_efficacy:.2%}\n")

print("Step 3: Compare the estimate with the true value")
print(f"Estimated Efficacy from Rates: {estimated_efficacy:.2%}")
print(f"True Per-Exposure Efficacy:      {true_efficacy:.2%}\n")

print("Conclusion: The two values are identical. Therefore, for an all-or-nothing vaccine,")
print("1 - Incidence Rate Ratio correctly estimates the per-exposure vaccine efficacy.")

<<<C>>>