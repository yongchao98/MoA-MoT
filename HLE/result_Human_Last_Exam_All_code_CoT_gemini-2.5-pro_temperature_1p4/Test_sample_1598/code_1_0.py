# The number of leaf vertices in the K_1,n graph.
n = 100

# We use a greedy approach to build the set of labels.
# We start with the first label x_1 = 1.
# Each subsequent label x_k is chosen to be the smallest integer that is not a sum
# of a subset of the preceding labels {x_1, ..., x_{k-1}}.
# This method generates labels that are powers of 2.
# x_1 = 1 = 2^0
# x_2 = 2 = 2^1
# x_3 = 4 = 2^2
# ...
# x_n = 2^(n-1)

# So, for n=100, the maximum label is 2^(100-1) = 2^99.

base = 2
exponent = n - 1

# Calculate the result
result = base ** exponent

# Print the final number as requested.
# The equation for the global labeling number k is k = base^exponent.
print(f"The global labeling number is {base}^{exponent}, which is:")
print(result)
