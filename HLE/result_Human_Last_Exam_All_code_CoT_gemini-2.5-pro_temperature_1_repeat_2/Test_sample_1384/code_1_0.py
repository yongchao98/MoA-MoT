import decimal

# We need high precision for this calculation. A precision of 50 decimal places
# is more than sufficient to ensure the accuracy of the first few digits.
decimal.getcontext().prec = 50

# The number we are evaluating is N = e^(-7^13).
# We express this in base 10 as N = 10^x, where x = -7^13 / ln(10).
# The equation for the digits is derived from N = 10^g * 10^K.

print("Calculating the first two non-zero digits of e^(-7^13).")
print("This is based on the scientific notation N = M * 10^K.")
print("-" * 20)

# Step 1: Calculate 7^13
power_val = decimal.Decimal(7) ** 13
print(f"The value of 7^13 is: {power_val}")

# Step 2: Calculate the exponent x = -7^13 / ln(10)
ln_10_val = decimal.Decimal(10).ln()
x = -power_val / ln_10_val
print(f"The full base-10 exponent x is: {x}")

# Step 3: Decompose x into its integer part K and fractional part g
# K is the exponent in scientific notation.
K = x.to_integral_value(rounding=decimal.ROUND_FLOOR)
# g is used to find the coefficient M.
g = x - K
print(f"The integer part of the exponent, K, is: {K}")
print(f"The fractional part of the exponent, g, is: {g}")


# Step 4: Calculate M = 10^g, which contains the non-zero digits.
M = decimal.Decimal(10) ** g
print(f"The coefficient M = 10^g is: {M}")
print("-" * 20)

# Step 5: Extract the first two digits from M.
# We convert M to a string for easy digit access. The first digit is before
# the decimal point, and the second is after.
m_str = f"{M:.{decimal.getcontext().prec-1}f}"
first_digit = m_str[0]
second_digit = m_str[2]  # Index 1 is the decimal point '.'

print(f"The first non-zero digit is the first digit of M: {first_digit}")
print(f"The second non-zero digit is the second digit of M: {second_digit}")

result = first_digit + second_digit
print(f"\nThus, the first two non-zero digits of e^(-7^13) are {result}.")
