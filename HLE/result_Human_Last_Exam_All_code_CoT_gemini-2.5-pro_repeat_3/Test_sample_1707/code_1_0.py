import math

# Define the given parameters as strings to handle the large exponents
alpha = "10^10000"
x0 = "10^-5000000"
y0 = "10^-5000000"

# The derived equation for T is:
# T = alpha / (x0/(1-x0) + y0/(1-y0))
# Since x0 = y0, this simplifies to T = alpha * (1 - x0) / (2 * x0)

# We print the final equation with the numbers substituted, as requested.
# The format clearly shows each number used in the calculation.
print("The final equation for T is:")
final_equation = f"T = ({alpha}) * (1 - {x0}) / (2 * {x0})"
print(final_equation)

# The equation can be simplified algebraically to find the exact value.
# T = (10^10000 * (1 - 10^-5000000)) / (2 * 10^-5000000)
# T = 0.5 * 10^(10000 - (-5000000)) * (1 - 10^-5000000)
# T = 0.5 * 10^5010000 * (1 - 10^-5000000)
# T = 0.5 * 10^5010000 - 0.5 * 10^(5010000 - 5000000)
# T = 0.5 * 10^5010000 - 0.5 * 10^10000

# We print the final, exact value of T in its symbolic form, as direct computation is infeasible.
final_value = "0.5 * 10^5010000 - 0.5 * 10^10000"
print("\nThe simplified exact value for T is:")
print(final_value)

# For clarity, this can also be written as:
final_value_alt = "5 * 10^5009999 - 5 * 10^9999"
print("\nAlternatively:")
print(final_value_alt)