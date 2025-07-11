import numpy as np

# Step 1: Define constants from the problem statement
# Half-life of Ba-140 (Parent) in days
T_half_Ba = 12.75
# Half-life of La-140 (Daughter) in days
T_half_La = 1.678
# Activity of La-140 at the first measurement (time T) in kBq/mL
A1 = 1.4
# Activity of La-140 at the second measurement (time T + 14 days) in kBq/mL
A2 = 2.1
# Time between measurements in days
dt = 14.0

# Step 2: Calculate the decay constants (lambda)
lambda_Ba = np.log(2) / T_half_Ba
lambda_La = np.log(2) / T_half_La

# Step 3: Solve for the time T
# The relationship derived from the Bateman equations is:
# exp((lambda_La - lambda_Ba) * T) = (A2 - A1 * exp(-dt * lambda_La)) / (A2 - A1 * exp(-dt * lambda_Ba))
# We will calculate the terms in this equation step by step.

# Calculate the decay factors over the 14-day period
decay_factor_La = np.exp(-dt * lambda_La)
decay_factor_Ba = np.exp(-dt * lambda_Ba)

# Calculate the numerator and denominator of the ratio
numerator = A2 - A1 * decay_factor_La
denominator = A2 - A1 * decay_factor_Ba

# Calculate the ratio
ratio = numerator / denominator

# Calculate the difference in decay constants
lambda_diff = lambda_La - lambda_Ba

# Finally, solve for T
T_days = np.log(ratio) / lambda_diff

# Step 4: Print the calculation and the result
print("--- Calculation of the time from irradiation to first analysis ---")
print(f"Half-life of Ba-140: {T_half_Ba} days")
print(f"Half-life of La-140: {T_half_La} days")
print(f"Decay constant of Ba-140 (λ_Ba): {lambda_Ba:.4f} day^-1")
print(f"Decay constant of La-140 (λ_La): {lambda_La:.4f} day^-1")
print("\nFirst La-140 activity (A1): 1.4 kBq/mL")
print("Second La-140 activity (A2): 2.1 kBq/mL")
print("Time between measurements (dt): 14 days")

print("\nThe time 'T' is calculated using the formula:")
print("T = ln( (A2 - A1 * exp(-λ_La * dt)) / (A2 - A1 * exp(-λ_Ba * dt)) ) / (λ_La - λ_Ba)")
print("\nPlugging in the values:")
print(f"T = ln( ({A2} - {A1} * exp(-{lambda_La:.4f} * {dt})) / ({A2} - {A1} * exp(-{lambda_Ba:.4f} * {dt})) ) / ({lambda_La:.4f} - {lambda_Ba:.4f})")
print(f"T = ln( ({A2} - {A1 * decay_factor_La:.4f}) / ({A2} - {A1 * decay_factor_Ba:.4f}) ) / ({lambda_diff:.4f})")
print(f"T = ln( {numerator:.4f} / {denominator:.4f} ) / {lambda_diff:.4f}")
print(f"T = ln({ratio:.4f}) / {lambda_diff:.4f}")
print(f"T = {np.log(ratio):.4f} / {lambda_diff:.4f}")
print(f"T = {T_days:.2f} days")
print("\nThe approximate time between sample irradiation and the first analysis is about 1.0 day.")

# Final answer in the required format
final_answer = round(T_days, 1)
print(f"\n<<< {final_answer} >>>")