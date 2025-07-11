import mpmath

# 1. Set the precision for the calculation.
# The number 10^100 has 101 digits. To find the remainder when divided by pi,
# we need to know pi with a precision greater than 101 digits.
# We set the precision (dps - decimal places) to 110, which is sufficient.
mpmath.mp.dps = 110

# 2. Define the number x = 10^100 as a high-precision number.
x = mpmath.mpf(10)**100

# 3. Use the periodicity of the tangent function, tan(x) = tan(x mod pi).
# We calculate y = 10^100 mod pi.
y = mpmath.fmod(x, mpmath.pi)

# 4. Calculate the tangent of the remainder y.
result = mpmath.tan(y)

# 5. Output the numbers involved in the final equation.
# To make the output readable, we format the high-precision numbers to 15 decimal places for printing.
y_str = mpmath.nstr(y, 15)
result_str = mpmath.nstr(result, 15)
print(f"tan(10^100) = tan(10^100 mod pi) = tan({y_str}...) = {result_str}...")

# 6. Extract and print the first three digits after the comma.
# The full result is a string like '-0.301551...'.
# We split the string representation at the decimal point and take the first 3 characters of the fractional part.
fractional_part = str(result).split('.')[1]
first_three_digits = fractional_part[:3]

print(f"\nThe first 3 digits after the comma of tan(10^100) are: {first_three_digits}")