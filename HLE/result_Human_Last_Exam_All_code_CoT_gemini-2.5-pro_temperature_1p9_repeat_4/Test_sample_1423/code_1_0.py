# The maximum number of distinct digits allowed.
k = 5

# The maximum length is given by the formula L_k = 2^k - 1.
# We calculate this value for k=5.
two_raised_to_k = 2**k
max_length = two_raised_to_k - 1

# Output the numbers in the final equation as requested.
print(f"Let k be the number of distinct digits. The maximum length of N is given by the formula 2^k - 1.")
print(f"For k = {k}, the equation is:")
print(f"2**{k} - 1 = {two_raised_to_k} - 1 = {max_length}")
print(f"\nThe maximum possible number of digits in N is {max_length}.")