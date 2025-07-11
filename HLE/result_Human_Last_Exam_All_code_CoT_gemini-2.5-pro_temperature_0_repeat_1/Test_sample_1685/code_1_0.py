import math

# Step 1: Define the equation for the coefficient C.
# From the dominant balance analysis, we have the equation C^3 = -6.

# Step 2: Solve for C.
# We need the real cube root of -6.
# C = (-6)^(1/3) = - (6)^(1/3)
C = -math.pow(6, 1/3)

# Step 3: Round the coefficient to two decimal places as requested.
C_rounded = round(C, 2)

# Step 4: The power of x, n, was found to be -1.
n = -1

# Step 5: Construct and print the final analytical expression.
# The expression is y(x) = C / x.
# The numbers in this equation are the coefficient C and the power n.
print(f"The approximate analytical expression for y(x) in the large x regime is of the form: y(x) = C * x^n")
print(f"The value for the power is n = {n}")
print(f"The value for the coefficient is C = {C_rounded:.2f}")
print("\nThe final expression is:")
print(f"y(x) = {C_rounded:.2f} / x")
