import sympy as sp

# Define the symbols used in the equation
K0, a, y, d, mu0, mu = sp.symbols('K_0 a y d mu_0 mu', real=True)
# Define i_x as a symbol for printing the vector direction
i_x = sp.Symbol('i_x')

# Based on the derivation, construct the terms of the final equation
# The constant factor in the formula
coefficient = mu0 / 2

# The term related to the current sheet amplitude
current_term = K0**2

# The term describing the spatial variation along the y-axis
y_variation_term = sp.sin(a * y)**2

# The denominator, which includes effects from the air gap and magnetic material
denominator_term = (sp.cosh(a * d) + (mu0 / mu) * sp.sinh(a * d))**2

# Combine the terms to form the magnitude of the force per area
force_magnitude = coefficient * current_term * y_variation_term / denominator_term

# The final expression is a vector in the x-direction
final_force_expression = force_magnitude * i_x

# Print the derived result in a structured way
print("The final expression for the force per unit area on the conductor is:")
sp.pprint(final_force_expression)

print("\nTo meet the requirement of outputting each number, let's look at the components:")
print(f"The numerical constant in the equation is: 1/2 (as in {sp.pretty(coefficient)})")
print("The other components are symbolic as shown in the full expression above.")
print("\nThis derived formula matches answer choice C.")
