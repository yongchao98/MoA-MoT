# This script determines and prints the expression for the relationship
# between 3-Hydroxypropionate ([B]) and PEP ([F]).

# The direct pathway is:
# 3-Hydroxypropionate -k2-> Malonyl-CoA -k3-> Acetyl-CoA -k4-> Pyruvate -k5-> PEP
# The concentration of the final product [F] is proportional to the initial
# reactant [B] and the rate constants of the direct path.

# Define the components of the equation
initial_reactant = "[B]"
final_product = "[F]"
# The numbers in the final equation are the subscripts of the rate constants
rate_constant_numbers = [2, 3, 4, 5]

# Build the expression string
expression_parts = [f"k{n}" for n in rate_constant_numbers]
full_expression = f"{final_product} ∝ {initial_reactant} * {' * '.join(expression_parts)}"

# Print the final equation
print("The expression that most accurately represents the relationship is:")
print(full_expression)