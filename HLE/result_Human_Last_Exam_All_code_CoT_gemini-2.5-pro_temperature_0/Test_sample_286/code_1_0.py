import math

# The problem reduces to the inequality 2*L >= C, where L = 1 - 10^{-M}
# and C is the supremum of the minimum possible range of partial sums.
# From discrepancy theory, we know C = 1.
# So, the inequality is 2 * (1 - 10^{-M}) >= 1.

# This simplifies to 1 - 10^{-M} >= 1/2.
# Which further simplifies to 1/2 >= 10^{-M}.
# Taking log base 10 on both sides, we get log10(1/2) >= -M.
# This is equivalent to -log10(2) >= -M, which means M >= log10(2).

# We now find the smallest positive integer M that satisfies this condition.

# The final inequality is M >= log_10(2)
val1 = 10
val2 = 2
print(f"The final inequality is M >= log_{val1}({val2})")

# Calculate the value of log_10(2)
log_val = math.log10(val2)
print(f"The value of log_{val1}({val2}) is approximately {log_val}")

# Find the smallest positive integer M
# Since M must be a positive integer, we take the ceiling of the log value.
# If the log value was 0 or negative, the smallest positive integer would still be 1.
if log_val <= 0:
    M = 1
else:
    M = math.ceil(log_val)

print(f"The smallest positive integer M satisfying the condition is {M}.")
