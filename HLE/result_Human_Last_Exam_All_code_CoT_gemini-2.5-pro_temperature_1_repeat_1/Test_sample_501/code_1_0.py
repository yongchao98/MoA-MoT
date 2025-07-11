import math

# Define the numerical constants present in the final derived formula.
# The derivation shows that the force is proportional to a coefficient,
# and inversely proportional to n and l raised to certain powers.
coefficient = 2
n_exponent = 2
l_exponent = 2

# Print the final force law using these constants.
# This code fulfills the requirement of outputting each number in the final equation
# by explicitly printing them as part of the formula string.
# Here, 'l' is used to represent the strut length \ell, and E(0) is the
# kinetic energy at zero extension.

print("The force law F(x) for a thermally isolated polymer chain at small extension x is:")
print(f"F(x) = - ({coefficient} * E(0) * x) / (n^{n_exponent} * l^{l_exponent})")
