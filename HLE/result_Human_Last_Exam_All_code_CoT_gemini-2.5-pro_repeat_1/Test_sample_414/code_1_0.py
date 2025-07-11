import numpy as np

# Half-lives of the nuclides in days
T_half_Ba = 12.75  # Half-life of Ba-140
T_half_La = 1.678  # Half-life of La-140

# Decay constants (lambda = ln(2) / T_half)
lambda_Ba = np.log(2) / T_half_Ba
lambda_La = np.log(2) / T_half_La

# The problem states that the activity increases from 1.4 to 2.1 kBq/mL.
# A direct simulation of the described experiment (pure Ba-140 decaying for 14 days)
# does not yield these results.
# A plausible interpretation is that the ratio of the total activity (Ba+La)
# to the parent activity (Ba) at the time of analysis was 2.1 / 1.4 = 1.5.
# This ratio can be used to determine the time (T) since the Ba-140 was pure
# (i.e., the cooling time after irradiation).
#
# The equation is: (A_Ba(T) + A_La(T)) / A_Ba(T) = 1.5
# This simplifies to: A_La(T) / A_Ba(T) = 0.5
#
# The ratio of activities A_La/A_Ba is given by:
# R(T) = (lambda_La / (lambda_La - lambda_Ba)) * (1 - exp(-(lambda_La - lambda_Ba) * T))
# We need to solve R(T) = 0.5 for T.

target_ratio = 0.5

# Rearranging the equation to solve for T:
# target_ratio = (lambda_La / (lambda_La - lambda_Ba)) * (1 - exp_term)
# exp_term = 1 - target_ratio * (lambda_La - lambda_Ba) / lambda_La
# - (lambda_La - lambda_Ba) * T = log(exp_term)
# T = -log(exp_term) / (lambda_La - lambda_Ba)

ratio_lambdas = lambda_La / (lambda_La - lambda_Ba)
exp_term = 1 - (target_ratio / ratio_lambdas)
T = -np.log(exp_term) / (lambda_La - lambda_Ba)

print("--- Calculation Steps ---")
print(f"Half-life of Ba-140: {T_half_Ba} days")
print(f"Half-life of La-140: {T_half_La} days")
print(f"Decay constant for Ba-140 (λ_Ba): {lambda_Ba:.4f} day^-1")
print(f"Decay constant for La-140 (λ_La): {lambda_La:.4f} day^-1")
print("\nInterpreting the data as a ratio of activities at the time of analysis:")
print("(Total Activity) / (Parent Activity) = 2.1 / 1.4 = 1.5")
print("This implies (Daughter Activity) / (Parent Activity) = 0.5")
print("\nSolving the equation for time T:")
print(f"(λ_La / (λ_La - λ_Ba)) * (1 - exp(-(λ_La - λ_Ba) * T)) = {target_ratio}")
print(f"{ratio_lambdas:.4f} * (1 - exp(-({lambda_La:.4f} - {lambda_Ba:.4f}) * T)) = {target_ratio}")
print(f"1 - exp(-{lambda_La - lambda_Ba:.4f} * T) = {target_ratio / ratio_lambdas:.4f}")
print(f"exp(-{lambda_La - lambda_Ba:.4f} * T) = {exp_term:.4f}")
print(f"-{lambda_La - lambda_Ba:.4f} * T = ln({exp_term:.4f}) = {np.log(exp_term):.4f}")
print(f"T = {-np.log(exp_term):.4f} / {lambda_La - lambda_Ba:.4f}")
print(f"T = {T:.2f} days")
print("\n--- Final Answer ---")
print(f"The approximate time between irradiation and the first analysis is {T:.2f} days.")