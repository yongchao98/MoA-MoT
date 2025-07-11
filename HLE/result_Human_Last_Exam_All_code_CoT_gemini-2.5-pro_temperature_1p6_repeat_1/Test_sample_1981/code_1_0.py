import math

# Step 1: Determine the value of X0 based on the problem's derivation.
# From the analysis of the solvability conditions and the integral equation,
# we find the relation: X0^(15/2) = 10^117.

# We solve for X0:
exponent_X0 = 117 * 2 / 15  # This simplifies to 15.6
X0 = 10**exponent_X0

# Step 2: Define the expression to be calculated.
# Expression: 10^30 * X0^2 - 10^30 * X0 + 10

# Step 3: Calculate each term of the expression.
# For printing clarity, we resolve the components of the expression.
term1_coeff = 10**30
term2_coeff = 10**30
constant_term = 10

term1_val = term1_coeff * (X0**2)
term2_val = term2_coeff * X0

# Step 4: Compute the final result.
result = term1_val - term2_val + constant_term

# Step 5: Print the final equation with each numerical part, as requested.
# Due to the large difference in magnitude, term2_val is negligible compared to term1_val
# in standard floating-point arithmetic, but we display it for completeness.
print("The final equation with numerical values is:")
print(f"{term1_val} - {term2_val} + {constant_term} = {result}")
