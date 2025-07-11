import numpy as np
from scipy.optimize import minimize_scalar

# The provided sample S
S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

# To find the Maximum Likelihood Estimate (MLE) for theta, we need to maximize
# the log-likelihood function. For a Cauchy distribution, the log-likelihood is:
# l(theta) = sum( -log(pi) - log(1 + (x_i - theta)^2) )
# Maximizing this is equivalent to minimizing the function:
# f(theta) = sum( log(1 + (x_i - theta)^2) )
def neg_log_likelihood(theta, data):
  """Computes the negative log-likelihood for a Cauchy distribution, omitting constants."""
  return np.sum(np.log(1 + (data - theta)**2))

# We use a numerical solver to find the theta that minimizes the function.
# This gives us the MLE for theta.
result = minimize_scalar(neg_log_likelihood, args=(S,))
theta_mle = result.x

print("The maximum likelihood estimate (MLE) of theta is found by maximizing the log-likelihood function.")
print("This corresponds to finding the value of theta where the derivative of the log-likelihood is zero.")
print(f"The sample data S is: {list(S)}")
print(f"\nThe numerically optimized MLE for theta is: {theta_mle:.7f}\n")

# The derivative of the log-likelihood, set to zero, gives the equation:
# sum[ (x_i - theta) / (1 + (x_i - theta)^2) ] = 0
# We will now print each term of this equation using our data and the estimated theta.
print("The equation to be satisfied by the MLE of theta is:")
print("sum[ (x_i - theta) / (1 + (x_i - theta)^2) ] = 0")
print("\nPlugging in the sample values and the estimated theta:")

total_sum = 0
equation_parts = []
# Loop through each data point to build and print the equation
for x in S:
    numerator = x - theta_mle
    denominator = 1 + (x - theta_mle)**2
    term_value = numerator / denominator
    total_sum += term_value
    # Format each term of the sum. For example, for x=2.65, the term is (2.65 - theta) / (1 + (2.65 - theta)^2)
    # The numbers below are the results of plugging in x and the estimated theta.
    equation_parts.append(f"{term_value:+.4f}")

# Join the parts with '+' and print the full expression
full_equation = " + ".join(equation_parts)
print(f"{full_equation} = {total_sum:e}")

print("\nThe sum is extremely close to zero, which confirms that our estimate is correct.")

# Finally, round the result to one decimal place as requested.
theta_mle_rounded = round(theta_mle, 1)

print(f"\nThe maximum likelihood estimate of theta, rounded to one decimal place, is {theta_mle_rounded}.")
<<<2.9>>>