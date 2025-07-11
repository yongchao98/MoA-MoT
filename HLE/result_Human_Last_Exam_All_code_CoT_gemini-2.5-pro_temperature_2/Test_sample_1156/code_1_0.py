import sympy
from sympy import sqrt, integrate, ln

# Define the symbolic variable
x = sympy.Symbol('x', real=True, positive=True)

# Define the unnormalized density function based on our derivation
# rho_unnormalized = 1 / (sqrt(x) * (1 + sqrt(x)))
# We will solve for the normalization constant C such that C * rho_unnormalized is the normalized density.
C = sympy.Symbol('C')
rho = C / (sqrt(x) * (1 + sqrt(x)))

# Calculate the integral of the unnormalized part to find the constant
integral_val = integrate(rho / C, (x, 0, 1))

# The normalization constant C is 1 / integral_val
C_val = 1 / integral_val
C_val_numeric = C_val.evalf()

# Create the normalized density function
rho_normalized = C_val / (sqrt(x) * (1 + sqrt(x)))

# Print the results
print("For the map T(x) = (1/sqrt(x) mod 1)^2:")
print("-" * 40)
print(f"The normalization constant is C = 1 / ({integral_val}) = {C_val}")
print(f"The numerical value of the constant is C ≈ {C_val_numeric}")
print("\nThe normalised density of the invariant measure is:")
sympy.pprint(rho_normalized, use_unicode=True)
print("\nFinal formula: rho(x) = (1 / (2 * ln(2))) / (sqrt(x) * (1 + sqrt(x)))")
# Output each number in the final equation
print(f"where 1 is one, 2 is two, ln(2) is approx {ln(2).evalf()},")
print(f"and the coefficient 1/(2*ln(2)) is approx {C_val_numeric}.")
