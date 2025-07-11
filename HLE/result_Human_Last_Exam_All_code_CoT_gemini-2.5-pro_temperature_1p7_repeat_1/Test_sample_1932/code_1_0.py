import math

# Define the symbols used in the final expression for the change in weight, Delta_W.
pi = "pi"
d = "d"
h = "h"
rho = "rho"
t = "t"
g = "g" # Note: g does not appear in the final correct expression.

# The derived expression for the change in weight is:
# Delta_W = (pi * d^2 * h^2 * rho) / (2 * t^2)

# According to the derived formula, we can break down the components.
# The expression contains the following variables raised to a power:
squared_vars = ["d", "h", "t"]
# The expression contains the following constants or variables:
linear_vars = ["rho", "pi"]
# The expression has a numerical coefficient in the denominator.
coefficient_denominator = 2

print("The change in weight, Delta_W, can be estimated by the expression corresponding to option C.")
print("The final equation is structured as a fraction.")
print(f"The numerator of the fraction consists of the product of several terms:")
print(f" - The constant '{pi}' (approximately {math.pi:.4f})")
print(f" - The variable '{d}' raised to the power of 2")
print(f" - The variable '{h}' raised to the power of 2")
print(f" - The variable '{rho}' raised to the power of 1")
print(f"The denominator of the fraction consists of the product of:")
print(f" - The number {coefficient_denominator}")
print(f" - The variable '{t}' raised to the power of 2")
print("\nFinal Expression: (pi * d**2 * h**2 * rho) / (2 * t**2)")