import sympy

# Define the symbol for the probability p
p = sympy.Symbol('p')

# Define the PDF of the Beta(2, 1) distribution
# f(p) = Gamma(2+1)/(Gamma(2)*Gamma(1)) * p^(2-1) * (1-p)^(1-1) = 2 * p
pdf = 2 * p

# Define the probability of hitting -1 for a given p
# H(p) = 1 if p <= 1/2
# H(p) = (1-p)/p if p > 1/2

# Calculate the first part of the integral from 0 to 1/2
integral1_val = sympy.integrate(1 * pdf, (p, 0, 1/2))

# Calculate the second part of the integral from 1/2 to 1
integral2_val = sympy.integrate(((1 - p) / p) * pdf, (p, 1/2, 1))

# The total probability is the sum of the two parts
total_prob = integral1_val + integral2_val

# Print the calculation steps and the final result
print("The probability is calculated as the sum of two integrals:")
print(f"1. Integral from 0 to 1/2 of (1 * 2p) dp = {integral1_val}")
print(f"2. Integral from 1/2 to 1 of ((1-p)/p * 2p) dp = {integral2_val}")
print(f"Total probability = {integral1_val} + {integral2_val} = {total_prob}")