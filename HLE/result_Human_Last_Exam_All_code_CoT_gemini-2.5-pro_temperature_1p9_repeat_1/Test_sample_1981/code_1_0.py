import math

# Calculate the exponent for X_0
# From the derivation, X_0^(15/2) = 10^117
# So, X_0 = 10^(117 * 2 / 15)
X0_exponent = 117 * 2 / 15

# Calculate X_0 and X_0^2
X0 = 10**X0_exponent
X0_sq = 10**(2 * X0_exponent)

# Define the coefficients of the final expression
coeff1 = 10**30
coeff2 = -1 * (10**30)
coeff3 = 10

# Calculate the final result
result = coeff1 * X0_sq + coeff2 * X0 + coeff3

# Print the components of the final equation and the result
# The numbers in the equation are the coefficients and the calculated values for X0 and X0^2
print(f"The equation to solve is: C1 * X0^2 + C2 * X0 + C3")
print(f"Where C1 = {coeff1}, C2 = {coeff2}, C3 = {coeff3}")
print(f"From the problem derivation, we found X0 = 10^{X0_exponent}")
print(f"X0 = {X0}")
print(f"X0^2 = {X0_sq}")
print("\nFinal Equation:")
# The instruction is to output each number in the final equation.
# This is ambiguous, it could mean showing the symbolic powers or the evaluated floats.
# We will show the floats.
print(f"{coeff1} * {X0_sq} - {abs(coeff2)} * {X0} + {coeff3} = {result}")

# The final result is a very large number, which can be printed in scientific notation.
print("\nFinal Answer:")
print(f"{result:.6e}")