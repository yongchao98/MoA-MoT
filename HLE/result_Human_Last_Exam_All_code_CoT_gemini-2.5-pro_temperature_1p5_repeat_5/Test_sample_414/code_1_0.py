import numpy as np
from scipy.optimize import fsolve

# --- Constants ---
# Half-life of Yttrium-90 in days
T_half_Y90 = 64.1 / 24.0  # 64.1 hours converted to days
# Decay constant (lambda) for Y-90 in day^-1
lambda_Y90 = np.log(2) / T_half_Y90

# Activity measurements
A1 = 1.4  # kBq/mL at time t
A2 = 2.1  # kBq/mL at time t + 14 days
ratio = A2 / A1

# --- Equation Setup ---
# We need to solve the equation:
# ratio = [1 - exp(-lambda * (t + 14))] / [1 - exp(-lambda * t)]
# for the unknown time 't'.

def equation_to_solve(t):
  """Defines the equation f(t) = 0 for the solver."""
  # t must be positive, handle t=0 case to avoid division by zero.
  if t <= 0:
    return -1 # Return a value that indicates it's not the root
  numerator = 1 - np.exp(-lambda_Y90 * (t + 14))
  denominator = 1 - np.exp(-lambda_Y90 * t)
  return numerator / denominator - ratio

# --- Solving the problem ---
# We make an initial guess for the solver. A few days seems reasonable.
initial_guess = 3.0
# Use fsolve to find the root of the equation
time_t = fsolve(equation_to_solve, initial_guess)[0]

# --- Output the results ---
print("Based on the analysis, the increase in activity is due to the ingrowth of Y-90 from its parent Sr-90.")
print("The time 't' between chemical separation and the first measurement can be found by solving the following equation:")
print(f"\n{A2}/{A1} = (1 - exp(-{lambda_Y90:.4f} * (t + 14))) / (1 - exp(-{lambda_Y90:.4f} * t))\n")
print(f"Which simplifies to:")
print(f"{ratio:.1f} = (1 - exp(-{lambda_Y90:.4f} * (t + 14))) / (1 - exp(-{lambda_Y90:.4f} * t))\n")

print(f"Solving for t gives the approximate time between separation and the first analysis.")
print(f"Assuming the time between irradiation and separation is negligible, this is our answer.\n")
print(f"Result: t = {time_t:.2f} days")

# Final answer format
print("\n<<<" + f"{time_t:.2f}" + ">>>")