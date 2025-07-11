# The number of distinct digits allowed is at most 5.
# To find the maximum possible number of digits in N, we should use the largest
# number of distinct digits, k.
k = 5

# The formula for the maximum length of a string where every subsequence
# has a unique character is 2^k - 1.
# We calculate this value for k=5.
max_length = 2**k - 1

# Print the final calculation, showing each number in the equation.
print(f"The maximum possible number of digits is given by the formula 2**k - 1.")
print(f"For k = {k}, the calculation is:")
print(f"2**{k} - 1 = {2**k} - 1 = {max_length}")