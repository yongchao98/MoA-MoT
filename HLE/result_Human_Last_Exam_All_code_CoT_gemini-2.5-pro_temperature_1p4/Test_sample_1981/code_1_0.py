import math

# Step 1: Define the value of X_0 based on the derivation.
# From the analysis, we found X_0^(15/2) = 10^117.
# Therefore, X_0 = 10^(117 * 2 / 15).
exponent_X0 = 117 * 2 / 15
X0 = 10**exponent_X0

# Step 2: Define the constants in the final expression.
C1 = 10**30
C2 = 10

# Step 3: Calculate the terms of the final expression and the result.
X0_squared = X0**2
term1 = C1 * X0_squared
term2 = C1 * X0
result = term1 - term2 + C2

# Step 4: Print the final equation with the computed numbers.
# The final equation is of the form C1 * X0^2 - C1 * X0 + C2
print("The final equation is: C1 * X0_squared - C1 * X0 + C2")
print("The numbers in the equation are:")
print(f"C1 = {C1}")
print(f"X0 = {X0}")
print(f"X0_squared = {X0_squared}")
print(f"C2 = {C2}")
print("\nFinal Calculation:")
# Using scientific notation for clarity
print(f"{C1:.1e} * {X0_squared:.4e} - {C1:.1e} * {X0:.4e} + {C2} = {result:.4e}")

# The final answer requested by the user is the result of the calculation
# Wrap the final numerical answer as requested.
print(f"\n<<< {result} >>>")