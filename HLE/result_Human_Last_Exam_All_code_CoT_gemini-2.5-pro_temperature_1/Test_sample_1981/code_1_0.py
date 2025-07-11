# Based on the derivation, X_0 is determined to be 1,
# assuming a typo correction in the problem's definition of beta.
X_0 = 1

# The expression to calculate is 10^30 * X_0^2 - 10^30 * X_0 + 10
# We define the coefficients and the constant term
coeff = 10**30
const = 10

# Calculate the result
result = coeff * X_0**2 - coeff * X_0 + const

# Print the equation with the final values
# Using integers for display where possible to show the cancellation clearly
print(f"The calculation is based on X_0 = {X_0}.")
print(f"The expression is 10**30 * X_0**2 - 10**30 * X_0 + 10.")
print(f"Substituting X_0 = {X_0}:")
# We show the equation with the value of X_0 plugged in
print(f"10**30 * {X_0}**2 - 10**30 * {X_0} + {const} = {result}")