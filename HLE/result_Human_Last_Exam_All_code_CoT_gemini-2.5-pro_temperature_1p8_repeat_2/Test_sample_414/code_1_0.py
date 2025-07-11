import numpy as np

# Step 1: Define constants
T_Ba = 12.75  # Half-life of Ba-140 in days
T_La = 1.678   # Half-life of La-140 in days
t_decay = 14.0 # Time between the two measurements in days

# Calculate decay constants (lambda) in days^-1
lambda_Ba = np.log(2) / T_Ba
lambda_La = np.log(2) / T_La

A1 = 1.4 # Activity at first measurement
A2 = 2.1 # Activity at second measurement
activity_ratio = A2 / A1

# Step 4: Solve for t1, the time from separation to the first measurement
# We need to solve the equation: activity_ratio = f(t1 + 14) / f(t1)
# where f(t) = exp(-lambda_Ba*t) - exp(-lambda_La*t)
# This is a transcendental equation, so we solve it numerically by searching.

def get_la_activity_ratio(t1):
    """Calculates the ratio of La-140 activity at t1+14 days to t1."""
    num = np.exp(-lambda_Ba * (t1 + t_decay)) - np.exp(-lambda_La * (t1 + t_decay))
    den = np.exp(-lambda_Ba * t1) - np.exp(-lambda_La * t1)
    if den == 0:
        return np.inf
    return num / den

# Numerical search for t1
# We search for t1 in a reasonable range (e.g., 0 to 5 days)
t1_solution = None
for t1_guess in np.arange(0.01, 5.0, 0.001):
    calculated_ratio = get_la_activity_ratio(t1_guess)
    if np.isclose(calculated_ratio, activity_ratio, atol=1e-4):
        t1_solution = t1_guess
        break
        
# Step 6: Calculate the cooling time (t_cool) based on the transient equilibrium assumption.
# We calculate the time to reach 99% of the transient equilibrium state.
f_equilibrium = 0.99
t_cool = -np.log(1 - f_equilibrium) / (lambda_La - lambda_Ba)

# Step 7: Calculate the final answer
total_time = t_cool + t1_solution

# Output the results and the final equation
print("Half-life of 140-Ba: {:.2f} days".format(T_Ba))
print("Half-life of 140-La: {:.3f} days".format(T_La))
print("Decay constant for 140-Ba (λ_Ba): {:.4f} days^-1".format(lambda_Ba))
print("Decay constant for 140-La (λ_La): {:.4f} days^-1".format(lambda_La))
print("-" * 30)
print("Time from separation to first measurement (t1): {:.2f} days".format(t1_solution))
print("Time to reach 99% transient equilibrium (t_cool): {:.2f} days".format(t_cool))
print("-" * 30)
print("The approximate time between irradiation and the first analysis is the sum of the cooling time and the processing time before the first measurement.")
print("Total Time = t_cool + t1")
print("Total Time = {:.2f} days + {:.2f} days = {:.2f} days".format(t_cool, t1_solution, total_time))
print("-" * 30)
# Round the final answer to one decimal place as it is an approximation
final_answer = round(total_time, 1)
print(f"Final Answer (approximate): {final_answer} days")
print(f"<<<{final_answer}>>>")