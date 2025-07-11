import math

# Define the variables as strings for formula representation
B = "B"
L = "log(1/delta)"

# Construct the formula string
# Note: In the context of the formula, 'log' refers to the natural logarithm.
# The formula combines the two identified asymptotic regimes.
asymptotic_formula = f"sqrt({B} * {L}) + {L} / log({L})"

print("The asymptotic value A(B, delta) such that d_{B,delta} = Theta(A(B, delta)) is given by:")
print(asymptotic_formula)
# The formula uses L = log(1/delta). Let's print it in terms of B and L as requested.
L_var = "L"
asymptotic_formula_with_L = f"sqrt({B} * {L_var}) + {L_var} / log({L_var})"
print("\nUsing the notation L = log(1/delta), the formula is:")
print(asymptotic_formula_with_L)