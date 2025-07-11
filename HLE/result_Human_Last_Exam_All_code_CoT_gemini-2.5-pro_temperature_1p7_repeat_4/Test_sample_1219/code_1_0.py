import mpmath

# 1. Set the precision for the calculation.
# The number 10^100 has 101 digits. To accurately find its remainder when
# divided by pi, our working precision must be greater than 101 digits.
# We will use 110 decimal places, which provides a safe margin.
mpmath.mp.dps = 110

# 2. Define the large number for the calculation.
x = mpmath.mpf(10)**100

# 3. Calculate the remainder and the final tangent value.
# The tangent function has a period of pi, so tan(x) = tan(x mod pi).
# We calculate this remainder using mpmath's high-precision operators.
rem = x % mpmath.pi
result = mpmath.tan(rem)

# 4. Print the final equation with the computed numbers.
# We format the high-precision numbers to a few decimal places for a clean display.
rem_str = mpmath.nstr(rem, 8)
result_str = mpmath.nstr(result, 8)
print("The calculation follows this equation:")
print(f"tan(10^100) = tan(10^100 mod pi) = tan({rem_str}...) = {result_str}...")

# 5. Extract and print the first 3 digits after the decimal point.
# To get the digits accurately, we convert the full-precision result to a string.
full_result_str = mpmath.nstr(result, mpmath.mp.dps)

# Find the location of the decimal point.
decimal_point_index = full_result_str.find('.')

# Extract the three characters immediately after the decimal point.
digits_after_comma = full_result_str[decimal_point_index + 1 : decimal_point_index + 4]

print(f"\nThe first 3 digits after the comma of tan(10^100) are: {digits_after_comma}")