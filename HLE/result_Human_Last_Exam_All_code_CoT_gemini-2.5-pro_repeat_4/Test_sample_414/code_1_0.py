import math

# --- Given Data ---
# Activity at first measurement in kBq/mL
A1 = 1.4
# Activity at second measurement in kBq/mL
A2 = 2.1
# Time between measurements in days
delta_t_days = 14
# Half-life of Yttrium-90 in hours
T_half_Y_hours = 64.1

# --- Calculations ---
# Convert time to consistent units (hours)
delta_t_hours = delta_t_days * 24

# Calculate the decay constant for Y-90 in hours^-1
lambda_Y_hours = math.log(2) / T_half_Y_hours

# We have the ratio of activities:
# A2 / A1 = [1 - exp(-lambda * (t + delta_t))] / [1 - exp(-lambda * t)]
# Let ratio = A2 / A1
# Let x = exp(-lambda * t)
# Let C = exp(-lambda * delta_t)
# The equation becomes: ratio = (1 - x*C) / (1 - x)
# Rearranging to solve for x: x = (ratio - 1) / (ratio - C)

ratio = A2 / A1
C = math.exp(-lambda_Y_hours * delta_t_hours)

# Solve for x
x = (ratio - 1) / (ratio - C)

# Now solve for t from x = exp(-lambda * t)
# -lambda * t = ln(x) => t = -ln(x) / lambda
t_hours = -math.log(x) / lambda_Y_hours

# Convert the final time to days
t_days = t_hours / 24

# --- Output Results ---
print("This script solves for the time 't' based on the ingrowth of a daughter nuclide.")
print("The governing equation is A2/A1 = [1 - exp(-lambda*(t+delta_t))] / [1 - exp(-lambda*t)]")
print("\n--- Constants and Parameters ---")
print(f"Half-life of Y-90 (T_1/2): {T_half_Y_hours} hours")
print(f"Decay constant of Y-90 (lambda): {lambda_Y_hours:.6f} hours^-1")
print(f"Activity A1: {A1} kBq/mL")
print(f"Activity A2: {A2} kBq/mL")
print(f"Time between measurements (delta_t): {delta_t_days} days ({delta_t_hours} hours)")

print("\n--- Final Equation with Solved Value 't' ---")
# This print statement shows how the solved 't' fits back into the ratio equation.
# It demonstrates the relationship between the measured activities and time.
verification_ratio = (1 - math.exp(-lambda_Y_hours * (t_hours + delta_t_hours))) / (1 - math.exp(-lambda_Y_hours * t_hours))
print(f"The solved equation with t = {t_hours:.2f} hours is:")
print(f"{A2} / {A1} = [1 - exp(-{lambda_Y_hours:.4f} * ({t_hours:.2f} + {delta_t_hours}))] / [1 - exp(-{lambda_Y_hours:.4f} * {t_hours:.2f})] = {verification_ratio:.2f}")


print("\n--- Answer ---")
print("The approximate time between the sample irradiation and the first analysis is:")
print(f"{t_days:.1f} days")

<<<4.2>>>