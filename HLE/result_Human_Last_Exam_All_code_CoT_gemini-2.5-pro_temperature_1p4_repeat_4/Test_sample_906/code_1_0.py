import math

# This script presents the final formula for the steady-state probability pi_0.
# The problem defines a birth-death process where:
# Birth rate from state i: lambda / (i+1)
# Death rate from state i: mu
# The traffic intensity is rho = lambda / mu.

# The derivation, as explained above, leads to the following result for pi_0.
# We will print the final equation and its components.

# Define the components of the formula as strings for clarity
base = "e"
exponent_operator = "^"
left_parenthesis = "("
negative_sign = "-"
rho_symbol = "rho"
right_parenthesis = ")"

print("The final derived formula for the steady-state probability pi_0 is:")

# We "output each number [or symbol] in the final equation" as requested
final_equation = f"pi_0 = {base}{exponent_operator}{left_parenthesis}{negative_sign}{rho_symbol}{right_parenthesis}"
print(final_equation)

print("\nWhere:")
print(f"  '{base}' is Euler's number, the base of the natural logarithm (approximately {math.e}).")
print(f"  '{rho_symbol}' is the traffic intensity, defined as rho = lambda / mu.")
