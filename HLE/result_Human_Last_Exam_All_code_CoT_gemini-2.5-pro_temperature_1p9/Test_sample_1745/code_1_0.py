import math

# Define the variables as strings for printing the formula.
# This code does not perform calculations but constructs the final expression string.
psi_y = "psi(y)"
z_1 = "z_1"
beta = "beta"
k = "k"
H = "H"
y = "y"

# Construct the final expression step-by-step
# Slip-dependent zeta potential at the bottom surface (j=1)
z_a1 = f"{z_1} * (1 + {beta} * {k})"

# Full expression for the EDL potential distribution
# Using the derived formula: psi(y) = z_a1 * sinh(k*(H-y)) / sinh(k*H)
# We will represent sinh(x) as "sinh(x)" in the string.
numerator = f"sinh({k} * ({H} - {y}))"
denominator = f"sinh({k} * {H})"
full_expression = f"{psi_y} = {z_a1} * {numerator} / {denominator}"

# Print the final expression
# The instruction was "Remember in the final code you still need to output each number in the final equation!"
# We interpret this as printing the symbolic form of the equation.
print("The expression for the Electrical double-layer potential distribution psi(y) is:")
print(f"psi(y) = z_1 * (1 + beta * k) * sinh(k * (H - y)) / sinh(k * H)")
