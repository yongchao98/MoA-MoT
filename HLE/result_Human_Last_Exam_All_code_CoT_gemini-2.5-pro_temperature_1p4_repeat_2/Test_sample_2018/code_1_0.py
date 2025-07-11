# The problem asks for the value of the computational factor used in the Carman-Kozeny
# equation in the original simulation-only implementation of the enthalpy-porosity method.
# This factor, A_mush, acts as a penalty term to suppress velocity in the solid phase.

# Based on historical accounts from the original authors, this value was chosen
# to be a large number before being refined in later benchmark studies.

# We define the numbers that make up the value in scientific notation.
base = 1
power = 10
exponent = 9

# Calculate the final value of the factor.
result = base * (power ** exponent)

# As requested, we print the numbers that form the final equation for the value.
print(f"The original computational factor was {base} * {power}^{exponent} = {int(result)}")