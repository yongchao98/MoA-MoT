import math

# Step 1: Define the constants for the Ba-140/La-140 decay chain.
# Half-life of the parent, Ba-140, in days.
T_p_days = 12.75
# Half-life of the daughter, La-140, in hours.
T_d_hours = 40.3
# Time between the two measurements in days.
delta_t_days = 14

# The ratio of the second measurement to the first.
activity_ratio = 2.1 / 1.4

# Step 2: Convert half-lives into decay constants (lambda = ln(2) / T).
# First, convert the daughter's half-life to days.
T_d_days = T_d_hours / 24
# Calculate decay constant for the parent (Ba-140).
lambda_p = math.log(2) / T_p_days
# Calculate decay constant for the daughter (La-140).
lambda_d = math.log(2) / T_d_days

# Step 3: Set up and solve the equation for the unknown time 't'.
# The ratio of daughter activities at t and t+delta_t is:
# Ratio = (exp(-lambda_p*(t+delta_t)) - exp(-lambda_d*(t+delta_t))) / (exp(-lambda_p*t) - exp(-lambda_d*t))
# Rearranging to solve for 't' gives:
# exp((lambda_d - lambda_p)*t) = (Ratio - exp(-lambda_d*delta_t)) / (Ratio - exp(-lambda_p*delta_t))

# Calculate the terms needed for the solution.
exp_lambda_p_delta_t = math.exp(-lambda_p * delta_t_days)
exp_lambda_d_delta_t = math.exp(-lambda_d * delta_t_days)

# Calculate the right-hand side (RHS) of the rearranged equation.
rhs_numerator = activity_ratio - exp_lambda_d_delta_t
rhs_denominator = activity_ratio - exp_lambda_p_delta_t
rhs = rhs_numerator / rhs_denominator

# Calculate the term on the left-hand side containing 't'.
ln_rhs = math.log(rhs)
lambda_diff = lambda_d - lambda_p

# Solve for 't', which is the time from irradiation to the first measurement.
t_initial_days = ln_rhs / lambda_diff

# Step 4: Print the final calculation and the result.
print("--- Calculation Steps ---")
print(f"Parent (Ba-140) half-life: {T_p_days} days")
print(f"Daughter (La-140) half-life: {T_d_hours} hours ({T_d_days:.3f} days)")
print(f"Parent decay constant (λ_p): {lambda_p:.5f} day⁻¹")
print(f"Daughter decay constant (λ_d): {lambda_d:.5f} day⁻¹")
print(f"Activity ratio (A2/A1): {activity_ratio}")
print("\nSolving the equation for t:")
print("exp((λ_d - λ_p) * t) = (A2/A1 - exp(-λ_d*14)) / (A2/A1 - exp(-λ_p*14))")
print(f"exp(({lambda_d:.5f} - {lambda_p:.5f}) * t) = ({activity_ratio} - {exp_lambda_d_delta_t:.5f}) / ({activity_ratio} - {exp_lambda_p_delta_t:.5f})")
print(f"exp({lambda_diff:.5f} * t) = {rhs_numerator:.5f} / {rhs_denominator:.5f}")
print(f"exp({lambda_diff:.5f} * t) = {rhs:.5f}")
print(f"{lambda_diff:.5f} * t = ln({rhs:.5f})")
print(f"{lambda_diff:.5f} * t = {ln_rhs:.5f}")
print(f"t = {ln_rhs:.5f} / {lambda_diff:.5f}")
print(f"t = {t_initial_days:.3f} days")
print("\n--- Final Answer ---")
print(f"The approximate time between irradiation and the first analysis is {t_initial_days:.2f} days.")
print(f"This is approximately {t_initial_days * 24:.1f} hours.")
<<<1.04>>>