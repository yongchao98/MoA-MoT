import sympy

# Define the symbols
pi, d, h, rho, t = sympy.symbols('pi d h rho t')

# Construct the expression for the change in weight based on the derivation
# Total mass of sand: M_s = rho * (pi * d**2 / 4) * h
# Acceleration of COM: a_COM = 2 * h / t**2
# Change in weight: delta_W = M_s * a_COM
delta_W = (rho * pi * d**2 * h / 4) * (2 * h / t**2)

# Simplify the expression
simplified_delta_W = sympy.simplify(delta_W)

# To match the format in the answer choices, we'll construct the string representation
# The simplified expression is (pi*d**2*h**2*rho)/(2*t**2)
# We need to print the formula using its components as requested
term_numerator = "pi * d**2 * h**2 * rho"
term_denominator = "2 * t**2"

print(f"The change in weight, Delta W, is given by the expression:")
print(f"Delta W = ({term_numerator}) / ({term_denominator})")

# In the format required for the final output, showing the numerator elements and denominator elements
# delta_W = pi * d^2 * h^2 * rho / (2 * t^2)
print("\nFinal Equation Form:")
print("pi * d**2 * h**2 * rho")
print("----------------------")
print("        2 * t**2")