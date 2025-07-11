from fractions import Fraction

# Step 1: Identify the exponent 'gamma' from the threshold N^(3/8).
gamma = Fraction(3, 8)

# Step 2: Choose the order of the moment 'p' for Chebyshev's inequality.
# Based on Bourgain's theorem for the L^p norm of the Schrodinger maximal function,
# the sharpest results are available for p=4.
p = 4

# Step 3: Determine the exponent 'beta' in the N-dependence of the integral bound.
# Bourgain's estimate, integral |M(x)|^4 dx <= C*N^(4*epsilon), implies that
# the effective exponent of N is 0 for calculating the main term of the bound.
beta = 0

# Step 4: Calculate alpha using the formula from Chebyshev's inequality.
# The inequality |X| * (N^gamma)^p <= integral |M(x)|^p dx leads to
# |X| <= N^(-p*gamma) * N^beta.
# Thus, the exponent alpha is beta - p * gamma.
alpha = beta - p * gamma

# Step 5: Print the final calculation and result.
# The problem asks to output each number in the final equation.
print("The exponent alpha is calculated using the formula: beta - p * gamma")
print("Substituting the values:")
final_equation = f"{beta} - {p} * {gamma.numerator}/{gamma.denominator} = {float(alpha)}"
print(final_equation)