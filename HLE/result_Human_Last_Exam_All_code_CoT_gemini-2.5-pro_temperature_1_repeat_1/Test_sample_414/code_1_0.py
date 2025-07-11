import math

# Half-lives of the nuclides
T_half_Ba140_days = 12.75  # Half-life of Ba-140 in days
T_half_La140_hours = 40.28 # Half-life of La-140 in hours

# Convert half-life of La-140 to days
T_half_La140_days = T_half_La140_hours / 24.0

# Calculate decay constants (lambda = ln(2) / T_half)
lambda_p = math.log(2) / T_half_Ba140_days  # Decay constant for parent Ba-140
lambda_d = math.log(2) / T_half_La140_days  # Decay constant for daughter La-140

# The problem is ambiguous as the cooling time before separation is unknown.
# A common interpretation for such problems is to calculate the characteristic time
# for the system, which is the time to reach maximum daughter activity (t_max)
# after separation. This represents the time to reach transient equilibrium.

# Calculate the time to reach maximum activity (t_max)
# Formula: t_max = ln(lambda_d / lambda_p) / (lambda_d - lambda_p)
t_max_numerator = math.log(lambda_d / lambda_p)
t_max_denominator = lambda_d - lambda_p
t_max_days = t_max_numerator / t_max_denominator

print("--- Radionuclide Properties ---")
print(f"Parent (Ba-140) half-life: {T_half_Ba140_days} days")
print(f"Daughter (La-140) half-life: {T_half_La140_hours} hours ({T_half_La140_days:.3f} days)")
print(f"Parent decay constant (λp): {lambda_p:.5f} day^-1")
print(f"Daughter decay constant (λd): {lambda_d:.5f} day^-1")
print("\n--- Calculation for Time to Maximum Activity (t_max) ---")
print(f"The time required for the daughter nuclide (La-140) to reach its maximum activity after separation is given by:")
print("t_max = ln(λd / λp) / (λd - λp)")
print(f"t_max = ln({lambda_d:.5f} / {lambda_p:.5f}) / ({lambda_d:.5f} - {lambda_p:.5f})")
print(f"t_max = ln({lambda_d / lambda_p:.3f}) / {t_max_denominator:.5f}")
print(f"t_max = {t_max_numerator:.5f} / {t_max_denominator:.5f}")
print(f"t_max = {t_max_days:.3f} days")
print("\n--- Conclusion ---")
print("Given the ambiguity in the problem statement (unknown cooling time), the most physically significant and calculable quantity is the time to reach maximum daughter activity.")
print("This characteristic time represents the point of transient equilibrium for the Ba-140/La-140 system.")
print(f"The approximate time is therefore {t_max_days:.1f} days.")
