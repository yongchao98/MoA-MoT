import math

# --- Given Data ---
# Basic reproduction number
R0 = 3.0
# Percentage of vaccinated people who still get the infection
infection_rate_in_vaccinated = 0.06

# --- Plan ---
# 1. Calculate Vaccine Effectiveness (E)
# The effectiveness is 1 minus the failure rate.
vaccine_effectiveness = 1 - infection_rate_in_vaccinated

# 2. Calculate the critical vaccination coverage (p_crit)
# The formula for herd immunity threshold is p_crit = (1 - 1/R0) / E
# This is the proportion of the population that needs to be effectively immune
# to prevent the epidemic spread (i.e., to get R_effective < 1).
critical_coverage_proportion = (1 - (1 / R0)) / vaccine_effectiveness

# Convert the proportion to a percentage for the final output
critical_coverage_percentage = critical_coverage_proportion * 100

# --- Output the Results ---
print("To prevent the spread of the virus, we need to calculate the theoretical threshold for vaccine coverage.")
print("This is based on the herd immunity formula: p_crit = (1 - 1/R₀) / E\n")

print(f"Given values:")
print(f"  - Basic Reproduction Number (R₀): {R0}")
print(f"  - Vaccine Effectiveness (E): 1 - {infection_rate_in_vaccinated} = {vaccine_effectiveness:.2f}\n")

print("Calculation of the required vaccine coverage:")
print(f"  p_crit = (1 - 1/{R0}) / {vaccine_effectiveness:.2f}")
print(f"  p_crit = ({1 - 1/R0:.4f}) / {vaccine_effectiveness:.2f}")
print(f"  p_crit = {critical_coverage_proportion:.4f}\n")

print(f"Therefore, the theoretical threshold of vaccine coverage we should aim for is {critical_coverage_percentage:.2f}%.")

# Final answer in the specified format at the end of the response
# The answer is the percentage value, rounded to two decimal places.
final_answer = round(critical_coverage_percentage, 2)