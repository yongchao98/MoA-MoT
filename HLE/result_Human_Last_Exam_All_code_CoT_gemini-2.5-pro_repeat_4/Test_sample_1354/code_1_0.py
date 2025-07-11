import numpy as np

# Problem parameters for the distributions
alpha = 3.0
beta = 2.0

# According to the properties of Gamma and Beta distributions,
# if a ~ Gamma(alpha, theta) and b ~ Gamma(beta, theta),
# then X = a / (a + b) follows a Beta(alpha, beta) distribution.
# The expectation of X is alpha / (alpha + beta).
mean_X = alpha / (alpha + beta)

# The random variable y_1 is defined as (a - b) / (a + b).
# This can be rewritten as 2 * a / (a + b) - 1 = 2*X - 1.
# Its expectation is E[y_1] = 2 * E[X] - 1.
E_y1 = 2 * mean_X - 1

# As derived in the explanation, the trace of the covariance matrix of v
# simplifies to 1 - (E[y_1])^2.
# Let's calculate this value.
trace_cov_v = 1 - E_y1**2

# --- Output the results ---
print("This script calculates the trace of the covariance matrix based on the analytical solution.")
print(f"The key intermediate value is E[(a-b)/(a+b)], which we call E_y1.")
print("\n--- Calculation Steps ---")
print(f"1. E[a/(a+b)] for a~Gamma({alpha},1), b~Gamma({beta},1) is {alpha}/({alpha}+{beta}) = {mean_X}")
print(f"2. E_y1 = 2 * E[a/(a+b)] - 1 = 2 * {mean_X} - 1 = {E_y1}")
print(f"3. The final equation is: Tr(Cov(v)) = 1 - (E_y1)^2")
print(f"   Plugging in the value of E_y1: 1 - ({E_y1})^2 = {trace_cov_v}")

print("\n--- Final Answer ---")
print(f"The trace of the covariance matrix is: {trace_cov_v}")
