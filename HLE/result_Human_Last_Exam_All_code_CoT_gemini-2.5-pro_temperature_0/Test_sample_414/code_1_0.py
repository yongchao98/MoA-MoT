import math

# Step 1: Define the physical constants for the Ba-140/La-140 system.
# Half-life of the parent, Ba-140 (in days)
T_p = 12.75
# Half-life of the daughter, La-140 (in days)
T_d = 40.22 / 24.0

# Calculate the decay constants (lambda) in units of days^-1
lambda_p = math.log(2) / T_p
lambda_d = math.log(2) / T_d

print("--- Physical Constants ---")
print(f"Parent (Ba-140) half-life (T_p): {T_p:.2f} days")
print(f"Daughter (La-140) half-life (T_d): {T_d:.2f} days")
print(f"Parent decay constant (λ_p): {lambda_p:.4f} days^-1")
print(f"Daughter decay constant (λ_d): {lambda_d:.4f} days^-1")
print("-" * 26)

# Step 2: Interpret the activity measurements from the problem.
# Assumed activity of the parent (Ba-140) at the time of separation (in kBq/mL)
A_Ba = 2.1
# Assumed activity of the daughter (La-140) at the time of separation (in kBq/mL)
A_La = 1.4

print("\n--- Interpreted Activities at Separation ---")
print(f"Activity of Ba-140 (A_Ba): {A_Ba} kBq/mL")
print(f"Activity of La-140 (A_La): {A_La} kBq/mL")
print("-" * 40)

# Step 3: Use the activity ratio to solve for the time since irradiation (t_wait).
# The ratio of activities, R, is a function of t_wait.
# R = A_La / A_Ba
R = A_La / A_Ba

# The equation for the ratio is: R = R_eq * (1 - exp(-(lambda_d - lambda_p) * t_wait))
# where R_eq is the ratio at transient equilibrium.
R_eq = lambda_d / (lambda_d - lambda_p)

# We rearrange the equation to solve for t_wait:
# t_wait = -ln(1 - R / R_eq) / (lambda_d - lambda_p)
try:
    # The term inside the logarithm must be positive
    log_term = 1 - (R / R_eq)
    if log_term <= 0:
        raise ValueError("The calculated activity ratio is not physically possible under this model.")

    t_wait = -math.log(log_term) / (lambda_d - lambda_p)

    print("\n--- Calculation Steps ---")
    print(f"Activity Ratio (R = A_La / A_Ba): {A_La} / {A_Ba} = {R:.4f}")
    print(f"Equilibrium Ratio (R_eq = λ_d / (λ_d - λ_p)): {lambda_d:.4f} / ({lambda_d:.4f} - {lambda_p:.4f}) = {R_eq:.4f}")
    print(f"Solving for t_wait in: {R:.4f} = {R_eq:.4f} * (1 - exp(-({lambda_d:.4f} - {lambda_p:.4f}) * t_wait))")
    print("-" * 25)

    print("\n--- Final Result ---")
    print(f"The approximate time between irradiation and analysis is: {t_wait:.2f} days.")
    print("-" * 20)

except ValueError as e:
    print(f"\nError: {e}")
