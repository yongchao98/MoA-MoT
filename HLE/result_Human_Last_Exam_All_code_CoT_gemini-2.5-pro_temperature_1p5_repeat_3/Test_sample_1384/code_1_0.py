import decimal

# Set the precision for decimal calculations. A precision of 50 digits is sufficient.
decimal.getcontext().prec = 50

# --- Step 1: Calculate 7^13 ---
power_val = 7**13
print(f"The equation we want to solve involves finding the first two digits of e^(-7^13).")
print(f"First, we calculate the exponent value:")
print(f"7^13 = {power_val}\n")

# --- Step 2: Calculate log10(e) ---
# log10(e) = 1 / ln(10)
ln_10 = decimal.Decimal(10).ln()
log10_e = decimal.Decimal(1) / ln_10
print(f"Next, we need the value of log10(e):")
print(f"log10(e) = {log10_e}\n")

# --- Step 3: Calculate log10(N) ---
# L = log10(N) = -7^13 * log10(e)
L = -decimal.Decimal(power_val) * log10_e
print(f"Now, we compute the base-10 logarithm of our number N = e^(-7^13):")
print(f"log10(N) = -7^13 * log10(e)")
print(f"log10(N) = {L}\n")

# --- Step 4: Isolate the fractional part ---
# k is the integer part (exponent in scientific notation)
# f is the fractional part, used to find the mantissa
k = L.to_integral_value(rounding=decimal.ROUND_FLOOR)
f = L - k
print(f"The number N can be written as 10^(k+f) = (10^f) * (10^k), where k is an integer and 0 <= f < 1.")
print(f"The integer part (exponent) is k = floor(log10(N)) = {k}")
print(f"The fractional part is f = log10(N) - k = {f}\n")

# --- Step 5: Compute the mantissa d = 10^f ---
d = decimal.Decimal(10)**f
print(f"The mantissa (significant digits) is d = 10^f:")
print(f"d = {d}\n")

# --- Step 6: Extract the first two digits ---
# The string representation of d will be "d.ddddd...".
# The first two non-zero digits are the digit before the decimal point and the first digit after.
d_str = str(d)
first_two_digits = d_str[0] + d_str[2]
print(f"The first two non-zero digits of e^(-7^13) are the first two digits of d.")
print(f"The first two digits are: {first_two_digits}")