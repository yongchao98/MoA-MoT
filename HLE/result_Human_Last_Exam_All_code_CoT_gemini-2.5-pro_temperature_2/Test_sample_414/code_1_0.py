import numpy as np
from scipy.optimize import fsolve

# Define the half-lives in days
T_p = 12.75  # Half-life of Ba-140 in days
T_d = 40.28 / 24.0  # Half-life of La-140 in days (40.28 hours)

# Calculate decay constants in days^-1
lambda_p = np.log(2) / T_p
lambda_d = np.log(2) / T_d

# Measured activities in kBq/mL
A0 = 1.4
A_t = 2.1

# The ratio of the activities
activity_ratio = A_t / A0

print("This script solves for the time 't' (in days) that satisfies the equation based on the Bateman equation for a freshly separated Ba-140 source.")
print("The equation describes the total activity A(t) relative to the initial parent activity A_p(0).")
print(f"We are solving for t in: A(t) / A_p(0) = {A_t} / {A0} = {activity_ratio:.3f}\n")
print(f"The equation to solve is:")
print(f"exp(-{lambda_p:.5f}*t) + ({lambda_d:.5f} / ({lambda_d:.5f} - {lambda_p:.5f})) * (exp(-{lambda_p:.5f}*t) - exp(-{lambda_d:.5f}*t)) = {activity_ratio:.3f}\n")


# Define the function to find the root of. We want to find t such that f(t) = 0.
def find_time(t):
    # This is the ratio A_total(t) / A_parent(0)
    term1 = np.exp(-lambda_p * t)
    k = lambda_d / (lambda_d - lambda_p)
    term2 = k * (np.exp(-lambda_p * t) - np.exp(-lambda_d * t))
    
    # We want to solve for when this expression equals the activity ratio
    return term1 + term2 - activity_ratio

# Use a numerical solver to find the time 't'. We start with an initial guess of 4 days.
# Based on preliminary analysis, the peak activity occurs around this time.
initial_guess = 4.0
solution_time = fsolve(find_time, initial_guess)[0]

print(f"The calculated time 't' that satisfies the equation is approximately {solution_time:.1f} days.")
print("This time represents the period between the initial separation/analysis and the point where the activity reached 2.1 kBq/mL, a 1.5-fold increase.")

# The final answer required by the user format
print("\nFinal Answer:")
print(f"{solution_time:.1f}")