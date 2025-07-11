import math

# Half-lives in days
T_p = 12.75  # Parent: Ba-140
T_d = 1.678  # Daughter: La-140

# Activity measurements in kBq/mL
A1 = 1.4
A2 = 2.1

# Time interval between measurements in days
delta_t = 14

# Calculate decay constants (lambda = ln(2) / T_half-life)
lambda_p = math.log(2) / T_p
lambda_d = math.log(2) / T_d

# Calculate the ratio of activities
R = A2 / A1

# The equation to solve for t1 is derived from the ratio of daughter activities:
# R = [exp(-lambda_p*(t1+delta_t)) - exp(-lambda_d*(t1+delta_t))] / [exp(-lambda_p*t1) - exp(-lambda_d*t1)]
# Rearranging gives:
# exp((lambda_d - lambda_p)*t1) = (R - exp(-delta_t*lambda_d)) / (R - exp(-delta_t*lambda_p))
# t1 = ln[ (R - exp(-delta_t*lambda_d)) / (R - exp(-delta_t*lambda_p)) ] / (lambda_d - lambda_p)

# Calculate terms for the solution
exp_term_p = math.exp(-delta_t * lambda_p)
exp_term_d = math.exp(-delta_t * lambda_d)

numerator = R - exp_term_d
denominator = R - exp_term_p

# Solve for t1, the time from separation to the first measurement
t1 = math.log(numerator / denominator) / (lambda_d - lambda_p)

print(f"Parent nuclide: Ba-140, Half-life = {T_p} days")
print(f"Daughter nuclide: La-140, Half-life = {T_d} days")
print(f"Decay constant of Ba-140 (λ_p): {lambda_p:.4f} day⁻¹")
print(f"Decay constant of La-140 (λ_d): {lambda_d:.4f} day⁻¹")
print("-" * 30)
print(f"First measurement (A1): {A1} kBq/mL")
print(f"Second measurement (A2) after {delta_t} days: {A2} kBq/mL")
print(f"Activity Ratio (R = A2/A1): {R:.4f}")
print("-" * 30)
print("The equation for the activity ratio is:")
print(f"{R:.1f} = [exp(-{lambda_p:.4f}*(t1 + {delta_t})) - exp(-{lambda_d:.4f}*(t1 + {delta_t}))] / [exp(-{lambda_p:.4f}*t1) - exp(-{lambda_d:.4f}*t1)]")
print("Solving this equation for t1 (the time from separation to the first analysis) gives:")
print(f"Time from separation to first analysis (t1) = {t1:.2f} days")
print("-" * 30)
print("Assuming the time between irradiation and chemical separation is negligible,")
print("the approximate time between irradiation and the first analysis is the value of t1.")
print(f"\nFinal Answer: The approximate time is {t1:.1f} days.")
