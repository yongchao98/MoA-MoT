import math

# Step 1: Identify parameters from the dominant term.
# The integrand is 1 / (epsilon + 9.0*x^5.0 + 5.0*x^6.0 + 9.0*x^8.0).
# For x near 0, the term 9.0*x^5.0 dominates the polynomial part.
# So, we have the form c*x^p.
c = 9.0
p = 5.0

# Step 2: Calculate the exponent for epsilon in the approximation.
# The scaling law for the integral gives the exponent E = (1-p)/p.
exponent = (1.0 - p) / p

# Step 3: Calculate the coefficient of the approximation.
# The coefficient C is given by the formula: c**(-1/p) * (pi/p) / sin(pi/p).
coefficient = (c**(-1.0 / p)) * (math.pi / p) / math.sin(math.pi / p)

# Step 4: Output the final analytical formula with the calculated numbers.
print("The analytical approximation for I(epsilon) for small epsilon is:")
print(f"I(epsilon) approx C * epsilon^E")
print("\nBased on the analysis, the constants are:")
print(f"C = {c}^(-1/{p}) * (pi/{p}) / sin(pi/{p})")
print(f"E = (1-{p})/{p}")
print("\nNumerically, these are:")
print(f"Coefficient C = {coefficient}")
print(f"Exponent E = {exponent}")
print("\nThus, the final approximate formula is:")
print(f"I(epsilon) \u2248 {coefficient} * epsilon^({exponent})")
