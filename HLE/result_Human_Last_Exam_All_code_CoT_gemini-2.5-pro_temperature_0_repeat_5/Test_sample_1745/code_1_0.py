import sympy

# Define the symbols used in the equation
# This allows us to represent the mathematical formula programmatically.
y, zeta_1, beta, k, H = sympy.symbols('y zeta_1 beta k H')

# The final derived expression for the Electrical double-layer (EDL) potential distribution.
# We construct this equation as a string for a clear and well-formatted output.
# This representation includes all variables and numbers (like 1 and 2 from H/2)
# as requested in the problem description.
final_equation = f"psi(y) = zeta_1*(1 + beta*k) * sinh(k*(H/2 - y)) / sinh(k*H)"

print("The final expression for the Electrical double-layer potential distribution is:")
print(final_equation)