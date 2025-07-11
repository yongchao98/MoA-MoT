import numpy as np

# The integral we need to compute is: integral from 0 to 1 of -3 * (2*e^(2x) + e^x) / (e^(2x) + e^x + 1) dx
# The analytical solution to this integral is -3 * [ln(e^(2x) + e^x + 1)] evaluated from x=0 to x=1.
# This results in the expression: -3 * (ln(e^2 + e + 1) - ln(3))

# Define the constant e
e = np.e

# Calculate the values of the terms inside the logarithms
upper_bound_term = e**2 + e + 1
lower_bound_term = 3

# Calculate the natural logarithms
ln_upper = np.log(upper_bound_term)
ln_lower = np.log(lower_bound_term)

# Calculate the final result of the definite integral
final_result = -3 * (ln_upper - ln_lower)

# Output the components of the final calculation
print("The final expression for the integral is: I = C * (ln(A) - ln(B))")
print(f"C = -3")
print(f"A = e^2 + e + 1 = {upper_bound_term}")
print(f"B = 3")
print("\nIntermediate calculations:")
print(f"ln(A) = ln({upper_bound_term}) = {ln_upper}")
print(f"ln(B) = ln({lower_bound_term}) = {ln_lower}")
print("\nFinal Result:")
print(f"I = -3 * ({ln_upper} - {ln_lower}) = {final_result}")
