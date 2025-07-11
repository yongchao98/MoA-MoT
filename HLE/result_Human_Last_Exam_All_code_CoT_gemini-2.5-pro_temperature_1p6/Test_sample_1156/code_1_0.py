import numpy as np

# Step 1: Define the structure of the normalized invariant density.
# The unnormalized density is rho_un(x) = 1 / (x + sqrt(x)).
# We need to find the normalization constant C such that C * rho_un(x) integrates to 1 over [0, 1].
# The integral is I = integral from 0 to 1 of 1/(x + sqrt(x)) dx.
# Using substitution u = sqrt(x), the integral becomes 2 * integral from 0 to 1 of 1/(u+1) du.
# This evaluates to 2 * [ln(u+1)] from 0 to 1 = 2 * (ln(2) - ln(1)) = 2 * ln(2) = ln(4).
# The normalization constant is therefore 1 / ln(4).

# Step 2: Calculate the numerical value of the constant ln(4).
constant = np.log(4)

# Step 3: Print the final equation for the normalized density.
# The format is rho(x) = 1 / (C * (x + sqrt(x))), where C is the calculated constant.
print("The normalised density of the invariant measure is given by the equation:")
print("rho(x) = 1 / (C * (x + sqrt(x)))")
print("\nwhere the constant C is ln(4).")
print("\nNumerically, C is approximately:")
print(constant)
print("\nSo the final equation is:")
print(f"rho(x) = 1 / ({constant} * (x + sqrt(x)))")

# Extract the numerical constant for the final answer block.
# The expression is rho(x) = 1/(ln(4)*(x+sqrt(x)))
final_answer_expression = f"1/({constant}*(x+sqrt(x)))"
