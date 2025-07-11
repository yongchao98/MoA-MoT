import math

# Given parameters from the problem description
alpha = 3.0
beta = 2.0

# The trace of the covariance matrix is derived to be (4 * alpha * beta) / (alpha + beta)^2
# This result is independent of d, theta, mu, Sigma, v1, and v2 under the given conditions.

# Calculate the numerator of the formula
numerator = 4 * alpha * beta

# Calculate the denominator of the formula
denominator = (alpha + beta)**2

# Compute the final result
trace_of_covariance = numerator / denominator

# Print the final equation with the substituted values and the result
print(f"The trace of the covariance matrix is calculated using the formula: (4 * alpha * beta) / (alpha + beta)^2")
print(f"Substituting the given values alpha = {alpha} and beta = {beta}:")
print(f"Trace = (4 * {alpha} * {beta}) / ({alpha} + {beta})^2 = {numerator} / {int(denominator)} = {trace_of_covariance}")
