# The number of distinct digits allowed.
# The problem states at most 5, and the length increases with the number of digits used,
# so we use k=5 for the maximum length.
k = 5

# The maximum length L(k) for a sequence with k distinct digits follows the recurrence L(k) = 2*L(k-1) + 1,
# with a base case L(1) = 1.
# This recurrence relation solves to L(k) = 2^k - 1.
# We calculate this value for k=5.
power_of_2 = 2**k
max_length = power_of_2 - 1

# Print the final result including the equation.
print("The maximum number of digits in N is determined by the formula L(k) = 2^k - 1, where k is the number of distinct digits.")
print(f"For at most 5 distinct digits, we use k = {k}.")
print("The calculation is:")
print(f"2^{k} - 1 = {power_of_2} - 1 = {max_length}")