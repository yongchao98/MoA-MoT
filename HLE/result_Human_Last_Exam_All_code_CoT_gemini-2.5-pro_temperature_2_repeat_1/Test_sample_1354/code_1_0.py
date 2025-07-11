from fractions import Fraction

# Given parameters from the problem description
alpha = 3
beta = 2

# Step 1: Calculate the expectation of V = a/(a+b), which follows a Beta(alpha, beta) distribution.
# The mean of a Beta(alpha, beta) distribution is alpha / (alpha + beta).
E_V_numerator = alpha
E_V_denominator = alpha + beta
E_V = Fraction(E_V_numerator, E_V_denominator)

# Step 2: Calculate the expectation of the first component of d, E[d_1].
# E[d_1] = E[(a-b)/(a+b)] = E[2*V - 1] = 2*E[V] - 1.
E_d1 = 2 * E_V - 1

# Step 3: Calculate the squared norm of the expectation of d, ||E[d]||^2.
# As shown in the derivation, E[d_i] = 0 for i > 1, so ||E[d]||^2 = (E[d_1])^2.
norm_E_d_squared = E_d1**2

# Step 4: Calculate the trace of the covariance matrix.
# Trace(Cov(v)) = Trace(Cov(d)) = 1 - ||E[d]||^2.
term1 = 1
term2 = norm_E_d_squared
trace_result = term1 - term2

# Print the step-by-step calculation
print("The trace of the covariance matrix is calculated as: Trace = 1 - ||E[d]||^2")
print(f"The parameters are: alpha = {alpha}, beta = {beta}")
print(f"The mean of V = a/(a+b) is E[V] = {E_V.numerator}/{E_V.denominator}")
print(f"The expectation of the first component of d is E[d_1] = 2 * E[V] - 1 = {E_d1.numerator}/{E_d1.denominator}")
print(f"The squared norm of the mean of d is ||E[d]||^2 = (E[d_1])^2 = ({E_d1.numerator}/{E_d1.denominator})^2 = {norm_E_d_squared.numerator}/{norm_E_d_squared.denominator}")

print("\nFinal equation and result:")
# Use the numbers from our calculation in the final equation string.
print(f"Trace = {term1} - {term2.numerator}/{term2.denominator} = {trace_result.numerator}/{trace_result.denominator}")
print(f"The numerical result is: {float(trace_result)}")
