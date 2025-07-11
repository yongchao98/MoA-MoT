import math

# Step 1: Define the parameters for the epsilon expansion.
# d is the spatial dimension of the system.
d = 3
# d_c is the upper critical dimension for the O(n) vector model.
d_c = 4
# n is the number of components of the order parameter.
# We will use n=1, which corresponds to the Ising model.
n = 1

# Step 2: Calculate the expansion parameter, epsilon.
epsilon = d_c - d

# Step 3: Calculate the specific heat scaling exponent alpha.
# The formula for alpha to the first order in epsilon is:
# alpha = (4 - n) / (2 * (n + 8)) * epsilon
numerator = 4 - n
denominator = 2 * (n + 8)
alpha = (numerator / denominator) * epsilon

# Step 4: Print the results and the calculation steps.
print("Calculating the specific heat scaling exponent 'alpha' using the epsilon expansion.")
print(f"Parameters used: d = {d} (spatial dimension), n = {n} (Ising model).")
print(f"The upper critical dimension is d_c = {d_c}.")
print(f"The expansion parameter is epsilon = d_c - d = {d_c} - {d} = {epsilon}.")
print("\nThe first-order formula for alpha is: alpha = (4 - n) / (2 * (n + 8)) * epsilon")

print("\nSubstituting the values into the formula:")
print(f"alpha = (4 - {n}) / (2 * ({n} + 8)) * {epsilon}")

# Show the step-by-step evaluation
print("\nEvaluating the expression:")
print(f"alpha = {numerator} / (2 * {n+8}) * {epsilon}")
print(f"alpha = {numerator} / {denominator} * {epsilon}")
print(f"alpha = {alpha}")

# Also print as a fraction
if isinstance(numerator, int) and isinstance(denominator, int):
    common_divisor = math.gcd(numerator, denominator)
    frac_num = numerator // common_divisor
    frac_den = denominator // common_divisor
    print(f"As a fraction, alpha = {frac_num}/{frac_den}")
