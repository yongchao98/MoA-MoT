import decimal

# Set the precision for decimal calculations. 50 digits is sufficient for this problem.
decimal.getcontext().prec = 50

print("To find the first two non-zero digits of e^(-7^13), we express the number in scientific notation d * 10^k.")
print("We need to find the first two digits of d.")
print("The calculation is based on the formula: d = 10^frac(-7^13 * log10(e))")
print("where frac(x) is the fractional part, x - floor(x).\n")

# Step 1: Calculate 7^13
val_7_13 = decimal.Decimal(7) ** 13
print(f"First, we calculate the exponent part:")
print(f"7^13 = {val_7_13}\n")

# Step 2: Calculate log10(e) with high precision. log10(e) = 1 / ln(10).
log10_e = decimal.Decimal(1) / decimal.Decimal(10).ln()
print(f"Next, we need the value of log10(e):")
print(f"log10(e) ≈ {log10_e}\n")

# Step 3: Compute the product -7^13 * log10(e)
log_val = -val_7_13 * log10_e
print(f"Now, we compute the full exponent for 10, which is -7^13 * log10(e):")
print(f"-{val_7_13} * {log10_e} ≈ {log_val}\n")

# Step 4: Find the fractional part of the result
# The floor of a number is the largest integer less than or equal to it.
integer_part = log_val.to_integral_value(rounding=decimal.ROUND_FLOOR)
fractional_part = log_val - integer_part
print(f"The fractional part of the result is:")
print(f"frac({log_val}) = {log_val} - ({integer_part}) ≈ {fractional_part}\n")

# Step 5: Calculate d = 10^(fractional part)
d = decimal.Decimal(10) ** fractional_part
print(f"Finally, we calculate d = 10^(fractional part):")
print(f"d = 10^{fractional_part} ≈ {d}\n")

# Step 6: Extract the first two non-zero digits from d
# d is between 1 and 10, so its string representation will be 'd.dddd...'.
d_str = str(d)
first_digit = d_str[0]
second_digit = d_str[2]  # The character at index 1 is the decimal point.
result = first_digit + second_digit
print(f"The first two non-zero digits of e^(-7^13) are the first two digits of d.")
print(f"The first two digits are {result}.")
