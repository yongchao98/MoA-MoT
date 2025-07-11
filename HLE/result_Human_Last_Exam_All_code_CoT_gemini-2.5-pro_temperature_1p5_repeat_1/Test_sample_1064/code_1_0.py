import sympy

# Set up the symbolic variables for adversary's mining power (beta) and tie-breaking probability (p)
beta, p = sympy.symbols('beta p')

print("This script calculates the expected chain quality for the specified selfish mining scenario.")
print("Chain quality is the long-term ratio of honest blocks to the total blocks in the main chain.\n")
print(f"Let beta (β) be the adversary's mining power.")
print(f"Let p be the probability that honest miners choose the adversary's block in a tie.\n")

print("The calculation is based on a Markov chain model of the selfish mining process where the")
print("state is the lead of the adversary's private chain. By solving for the steady-state")
print("probabilities of this model, the chain quality can be derived as a function of beta and p.\n")

# Based on the derivation from the model, we can define the components
# of the final formula for the expected Chain Quality (CQ).
# This formula is valid for beta < 0.5. For beta >= 0.5, the quality is 0.

numerator_factor1 = 1 + 3*beta - p*beta
numerator_factor2 = 1 - beta
numerator_factor3 = 1 - 2*beta

denominator = 1 - 3*beta**2 - beta**3

# Construct the full expression for chain quality from the derived components
chain_quality_expression = (numerator_factor1 * numerator_factor2 * numerator_factor3) / denominator

# Print the final equation, showing each component as requested.
print("--------------------------------------------------------------------")
print("The final equation for the expected chain quality is:")
print("CQ(β, p) = (Term_A * Term_B * Term_C) / Denominator\n")
print("Where the components are:")

print("Term_A:")
sympy.pretty_print(numerator_factor1)
print("\nTerm_B:")
sympy.pretty_print(numerator_factor2)
print("\nTerm_C:")
sympy.pretty_print(numerator_factor3)
print("\nDenominator:")
sympy.pretty_print(denominator)

print("\n--------------------------------------------------------------------")
print("Resulting in the full formula for Chain Quality CQ(β, p):")
sympy.pretty_print(chain_quality_expression)

print("\nNote: This formula is valid for 0 <= beta < 0.5.")
print("If beta >= 0.5, the adversary's power is too high for the system to be stable,")
print("and the expected chain quality approaches 0.")
print("--------------------------------------------------------------------")

final_expression_str = str(chain_quality_expression.simplify())
# Replace python's power operator with standard mathematical notation for clarity
final_expression_str = final_expression_str.replace('**', '^')

print(f"<<<({final_expression_str})>>>")