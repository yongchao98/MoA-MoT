# The basic reproduction number of the virus.
R0 = 3.0

# The percentage of vaccinated people who still get infected.
breakthrough_rate_percent = 6.0
breakthrough_rate = breakthrough_rate_percent / 100

# Calculate the vaccine's effectiveness (E) at preventing infection.
# E = 1 - (rate of infection in vaccinated individuals)
vaccine_effectiveness = 1 - breakthrough_rate

# The formula to calculate the critical vaccination coverage (p_c) is:
# p_c = (1 - 1/R0) / E
# This formula determines the proportion of the population that needs to be
# vaccinated to bring the effective reproduction number below 1.

# Calculate the critical vaccination coverage.
critical_coverage = (1 - 1/R0) / vaccine_effectiveness

# Convert the result to a percentage for readability.
critical_coverage_percent = critical_coverage * 100

# --- Output the results ---
print("--- Calculating Required Vaccine Coverage ---")
print(f"Basic Reproduction Number (R0): {R0}")
print(f"Vaccine Effectiveness (E): 1 - {breakthrough_rate} = {vaccine_effectiveness:.2f} or {vaccine_effectiveness*100:.0f}%")
print("\nThe formula for the critical vaccination coverage (p_c) is: p_c = (1 - 1/R0) / E")
print("\nPlugging in the numbers:")
print(f"p_c = (1 - 1/{R0}) / {vaccine_effectiveness:.2f}")
print(f"p_c = ({1 - 1/R0:.3f}) / {vaccine_effectiveness:.2f}")
print(f"p_c = {critical_coverage:.4f}")

print("\n--- Conclusion ---")
print(f"To prevent the spread of the virus, the theoretical threshold of vaccine coverage you should reach is {critical_coverage_percent:.1f}%.")
