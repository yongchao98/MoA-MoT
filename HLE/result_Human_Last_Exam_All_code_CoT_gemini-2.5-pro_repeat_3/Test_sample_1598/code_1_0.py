# The problem of finding the global labeling number of K_1,100 reduces to finding
# the minimum possible value for the largest label in a set of 100 positive integers
# {l_1, ..., l_100}, where each label l_i is not equal to any sum of a subset of the other 99 labels.

# As derived in the explanation, the optimal strategy to minimize the maximum label is to
# choose labels that are powers of 2. For 100 labels, the optimal set is
# {2**0, 2**1, 2**2, ..., 2**99}.

# The global labeling number is the maximum value in this set.
# The graph has n=100 leaf vertices.
n = 100

# The largest label will be 2**(n-1).
base = 2
exponent = n - 1

# Calculate the result.
result = base ** exponent

# Print the final equation with the numbers involved.
print(f"The global labeling number for K(1,100) is given by the equation:")
print(f"{base}**{exponent} = {result}")
