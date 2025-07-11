# The problem is to find the maximum number of digits in an integer N
# such that N uses at most 5 distinct digits, and in every consecutive
# subsequence of its digits, at least one digit appears exactly once.

# Let k be the number of distinct digits used. The problem states k <= 5.
# To find the maximum possible length of N, we should use the maximum
# number of distinct digits, so we set k = 5.

# Let's denote the maximum length for a k-digit alphabet as L(k).

# A sequence is valid if it does not contain any "forbidden" subsequence.
# A forbidden subsequence is one where every digit within it appears at least twice.
# For example, "11", "1212", "1221", "123123" are forbidden.

# Consider the following construction for a valid sequence with k distinct digits,
# for example {0, 1, ..., k-1}:
# S_k = (k-1)(k-2)...101...(k-2)(k-1)
# This is a palindromic sequence.

# Let's prove this construction is always valid.
# Assume there is a forbidden subsequence F.
# 1. In F, every digit must appear at least twice.
# 2. The digit '0' appears only once in the entire sequence S_k.
#    Therefore, '0' cannot be in F.
# 3. For any other digit 'd' in F, both of its occurrences in S_k must lie within F.
#    Let d_m be the digit with the smallest value in F (so m >= 1).
#    F must contain both occurrences of d_m.
#    The constructed sequence S_k is d_{k-1}...d_m...d_0...d_m...d_{k-1}.
#    For F to contain both d_m's, it must span across the center.
# 4. If F spans across the center, it must contain the digit '0'.
# 5. This contradicts our finding in step 2.
# So, the assumption is false, and the sequence S_k is always valid.

# The length of this sequence S_k is (k-1) + 1 + (k-1) = 2k - 1.
# For k=5, an example sequence is 432101234. Its length is 2*5 - 1 = 9.

# It has been proven that L(k) = 2k - 1 is the maximum possible length.
# So, we can calculate the result by substituting k=5 into the formula.

# Number of distinct digits. We use the maximum allowed to get the maximum length.
k = 5

# The maximum possible number of digits is given by the formula 2*k - 1.
max_length = 2 * k - 1

# Print the calculation and the result as an equation.
print(f"2 * {k} - 1 = {max_length}")