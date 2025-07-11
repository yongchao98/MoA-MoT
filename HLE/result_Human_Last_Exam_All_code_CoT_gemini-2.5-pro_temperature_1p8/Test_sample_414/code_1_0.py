import numpy as np
from scipy.optimize import fsolve

# Step 1 & 2: Define constants for the Ba-140/La-140 pair
T_half_parent_days = 12.75  # Half-life of Ba-140 in days
T_half_daughter_days = 40.22 / 24.0  # Half-life of La-140 in days (converted from hours)

# Calculate decay constants (lambda) in units of day^-1
lambda_p = np.log(2) / T_half_parent_days
lambda_d = np.log(2) / T_half_daughter_days

# Measured activities and time interval
A1 = 1.4  # kBq/mL
A2 = 2.1  # kBq/mL
ratio_A2_A1 = A2 / A1
delta_t_days = 14.0

# Step 3 & 4: Set up the equation to solve for t_prime
# A(t) is proportional to (exp(-lambda_p*t) - exp(-lambda_d*t))
# We need to solve ratio = A(t_prime + delta_t) / A(t_prime) for t_prime
def equation_to_solve(t_prime):
    """
    This function represents the equation ratio - f(t+14)/f(t) = 0.
    We are looking for the root, which is t_prime.
    t_prime must be a single float, not an array.
    """
    numerator = np.exp(-lambda_p * (t_prime + delta_t_days)) - np.exp(-lambda_d * (t_prime + delta_t_days))
    denominator = np.exp(-lambda_p * t_prime) - np.exp(-lambda_d * t_prime)
    
    # Avoid division by zero if t_prime is 0
    if np.abs(denominator) < 1e-9:
        return ratio_A2_A1 - 1e9 # Return a large number if t_prime is close to 0
        
    return ratio_A2_A1 - (numerator / denominator)

# Initial guess for the solver. Must be a small positive number.
initial_guess = 1.0
t_prime_solution = fsolve(equation_to_solve, initial_guess)
t_prime = t_prime_solution[0]

# Step 5: Assume cooling time is approximately one half-life of the parent
T_cooling = T_half_parent_days

# Step 6: Calculate the total approximate time
total_time = T_cooling + t_prime

# Output the results step-by-step
print(f"1. Identified radioactive pair: Ba-140 (Parent, T_1/2 = {T_half_parent_days:.2f} days) -> La-140 (Daughter, T_1/2 = {T_half_daughter_days:.2f} days).")
print(f"2. Calculated time from chemical separation to the first analysis (t'): {t_prime:.2f} days.")
print(f"3. Assumed cooling time from irradiation to separation (T_c) is one parent half-life: {T_cooling:.2f} days.")
print(f"4. The approximate time between irradiation and the first analysis is the sum of these two periods.")
print(f"   Total Time = T_c + t'")
print(f"   Total Time = {T_cooling:.2f} days + {t_prime:.2f} days = {total_time:.2f} days.")
print(f"\nThe final approximate time is {total_time:.1f} days.")
# The problem contains a subtle paradox in its numbers, but this solution follows the most physically plausible model and assumptions.
# Given the request for an approximate time, we provide a rounded value.
final_answer = round(total_time)
print(f"\nThus, the final approximate answer is {final_answer} days.")

# Return format as requested
# <<<14>>>