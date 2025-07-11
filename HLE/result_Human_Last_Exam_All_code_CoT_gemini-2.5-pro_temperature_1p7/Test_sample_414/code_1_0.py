import numpy as np
from scipy.optimize import fsolve

# Step 1: Define half-lives and calculate decay constants
T_Ba = 12.752  # Half-life of Ba-140 in days
T_La = 1.678   # Half-life of La-140 in days

lambda_Ba = np.log(2) / T_Ba
lambda_La = np.log(2) / T_La

# The equation to solve comes from the ratio of La-140 activities at two different times.
# A_La(t1 + 14) / A_La(t1) = 1.5
# Let the function D(t) be the time-dependent part of the La-140 activity equation:
# D(t) = exp(-lambda_Ba * t) - exp(-lambda_La * t)
# We need to solve: D(t + 14) / D(t) = 1.5, which is equivalent to D(t + 14) - 1.5 * D(t) = 0.

def equation_to_solve(t):
    """
    Represents the equation D(t + 14) - 1.5 * D(t) = 0, where t is the time from
    separation to the first measurement.
    """
    # The time t must be positive.
    if t <= 0:
        return np.inf
        
    d_t1 = np.exp(-lambda_Ba * t) - np.exp(-lambda_La * t)
    d_t2 = np.exp(-lambda_Ba * (t + 14)) - np.exp(-lambda_La * (t + 14))
    
    return d_t2 - 1.5 * d_t1

# Step 2: Numerically solve the equation for t.
# We need an initial guess for the solver. Based on the physics, the time should be
# on the order of a few days. Let's start with a guess of 1 day.
initial_guess = 1.0
time_to_first_analysis = fsolve(equation_to_solve, initial_guess)[0]

# Step 3: Print the final answer and the equation.
print("The approximate time between separation and the first analysis is:")
print(f"{time_to_first_analysis:.2f} days")
print("\nThis result is obtained by solving the following equation for t:")
print("(exp(-{:.3f}*(t+14)) - exp(-{:.3f}*(t+14))) / (exp(-{:.3f}*t) - exp(-{:.3f}*t)) = 1.5".format(lambda_Ba, lambda_La, lambda_Ba, lambda_La))
print(f"where t = {time_to_first_analysis:.2f} days.")
