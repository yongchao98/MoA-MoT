import decimal

# Set the precision for the calculation. 50 digits is sufficient.
decimal.getcontext().prec = 50

# Let N = e^(-7^13). We want to find its first two non-zero digits.
# We express N in scientific notation as N = K * 10^E, where 1 <= K < 10.
# The leading digits of N are the digits of K.
# We find K by computing the base-10 logarithm of N.
# log10(N) = -7^13 * log10(e) = -7^13 / ln(10)
log10_N = -(decimal.Decimal(7)**13) / decimal.Decimal(10).ln()

# The logarithm log10(N) can be split into an integer part E and a fractional part M:
# log10(N) = E + M, where E = floor(log10(N)) and 0 <= M < 1.
# The mantissa of N is then K = 10^M.
M = log10_N - log10_N.to_integral_value(rounding=decimal.ROUND_FLOOR)
K = decimal.Decimal(10)**M

# The final equation to find the mantissa is K = 10^M.
# The values of M and K are:
print(f"The fractional part of the logarithm (M) is: {M}")
print(f"The mantissa (K = 10^M) is: {K}")

# The first two non-zero digits of e^(-7^13) are the first two digits of K.
k_str = str(K)
first_two_digits = k_str[0] + k_str[2] # k_str[1] is the decimal point

print(f"The first two non-zero digits of e^(-7^13) are: {first_two_digits}")