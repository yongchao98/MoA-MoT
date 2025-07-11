import math

# Step 1: Calculate X_0 based on the derived equation X_0^(15/2) = 10^117
exponent_X0 = 117 * (2.0 / 15.0)
X0 = 10**exponent_X0

# Step 2: Define the expression to be calculated
# Expression: 10^30 * X_0^2 - 10^30 * X_0 + 10

# Step 3: Calculate each term of the expression
term1_coeff = 10**30
term2_coeff = 10**30
constant = 10

X0_squared = X0**2

term1 = term1_coeff * X0_squared
term2 = term2_coeff * X0

result = term1 - term2 + constant

# Step 4: Output the numbers in the final equation and the result
print("The value of X_0 is calculated as 10 raised to the power of 15.6:")
print(f"X_0 = {X0}")
print("\nThe expression to calculate is:")
print(f"{term1_coeff} * X_0^2 - {term2_coeff} * X_0 + {constant}")
print("\nSubstituting the calculated values into the expression:")
print(f"{term1} - {term2} + {constant} = {result}")
