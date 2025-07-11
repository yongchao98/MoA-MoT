import math

# This script calculates and prints the symbolic expression for the change
# in mutual inductance per unit length (Δm) between the two circuits
# when the concentrators are added.

# --- Symbolic representation of parameters ---
# These are the variables from the problem description.
# We use Unicode for a clear representation of Greek letters and subscripts.
mu_0 = "μ₀"  # Permeability of free space
pi = "π"    # The mathematical constant pi
h = "h"      # Separation between wires in a single circuit
d = "d"      # Distance between the centers of the two circuits
R1 = "R₁"    # Inner radius of the concentrator shells
R2 = "R₂"    # Outer radius of the concentrator shells

# --- The final derived expression ---
# The change in mutual inductance per unit length (Δm) is derived from:
# Δm = m₁ * (enhancement_factor - 1), where m₁ is the inductance for bare wires.
# For d >> h, the mutual inductance per unit length is m₁ ≈ - (μ₀ * h²) / (2 * π * d²).
# The total enhancement factor is the product of the source enhancement (R₂/R₁) 
# and the receiver enhancement (R₂/R₁), which gives (R₂/R₁)².

# We will now construct and print the equation string. The prompt requests to "output each number
# in the final equation", which we interpret as displaying all the symbols and constants
# that form the expression.
print("The expression for the change in mutual inductance per unit length (Δm) is:")

# Creating a formatted string for the final equation
final_equation = f"Δm = - ( {mu_0} * {h}² ) / ( 2 * {pi} * {d}² ) * [ ( {R2} / {R1} )² - 1 ]"

print(final_equation)
