import math

# The problem is to find the supremum of X, which we derived to be
# X_sup = 24 / (15 + 16 * pi^2)

# Define the numbers in the final equation
numerator = 24.0
denominator_term1 = 15.0
denominator_term2_coeff = 16.0

# Calculate pi squared
pi_squared = math.pi**2

# Calculate the denominator
denominator_value = denominator_term1 + denominator_term2_coeff * pi_squared

# Calculate the final value
supremum_X = numerator / denominator_value

print("The expression for the supremum of X is: 24 / (15 + 16 * \u03c0^2)")
print("\n--- Calculation Breakdown ---")
print(f"Numerator: {numerator}")
print(f"Value of \u03c0^2 (pi squared): {pi_squared}")
print(f"Denominator is {denominator_term1} + {denominator_term2_coeff} * {pi_squared}")
print(f"This evaluates to: {denominator_value}")
print("\n--- Final Answer ---")
print(f"Supremum of X = {numerator} / {denominator_value}")
print(f"Supremum of X = {supremum_X}")
