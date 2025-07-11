import decimal

# 1. Set the precision for decimal calculations.
# The number 10^100 has 101 digits. The precision must be higher than that
# to get an accurate remainder. We'll use 120.
decimal.getcontext().prec = 120

# 2. Use a high-precision string for Pi. This value has 125 decimal places.
PI_STR = "3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647"
pi = decimal.Decimal(PI_STR)

# 3. Define the number x = 10^100 as a Decimal object.
x = decimal.Decimal(10) ** 100

# 4. Calculate the remainder of x divided by pi.
# This is the core of the problem, using the periodicity of tan.
rem = x % pi

# 5. Calculate the tangent of the remainder.
tan_val = rem.tan()

# The final equation and its components:
print(f"tan(10^100) is equivalent to tan(10^100 mod π)")
print(f"The remainder (10^100 mod π) is: {rem}")
print(f"The final value of tan(10^100) is: {tan_val}")

# 6. Extract and print the first 3 digits after the decimal point.
# Convert the result to a string to easily access the digits.
tan_str = str(tan_val)
# Find the decimal point and slice the next 3 characters.
decimal_point_index = tan_str.find('.')
first_three_digits = tan_str[decimal_point_index + 1 : decimal_point_index + 4]

print(f"\nThe first 3 digits after the comma of tan(10^100) are: {first_three_digits}")