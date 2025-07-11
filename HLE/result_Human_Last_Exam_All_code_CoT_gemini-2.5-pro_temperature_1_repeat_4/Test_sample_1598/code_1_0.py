# The graph is K_1,n where n = 100.
n = 100

# As derived in the explanation, the problem reduces to finding a set of 100 positive
# integer labels {x_1, ..., x_100} such that no label is a sum of any subset
# of the other 99 labels. To find the global labeling number, we need to find
# such a set where the maximum label is minimized.

# A greedy construction where we choose the smallest possible integer at each step
# yields the sequence x_k = 2^(k-1) for k = 1, 2, ..., n.
# This is because the sums of subsets of {1, 2, 4, ..., 2^(k-2)} produce
# all integers from 1 to 2^(k-1) - 1. Therefore, the smallest integer
# not in this set of sums (and greater than the previous term 2^(k-2)) is 2^(k-1).

# The global labeling number is the maximum label in this set, which is the n-th term.
# The equation for the n-th term is x_n = 2^(n-1).

base = 2
exponent = n - 1

# Calculate the result
global_labeling_number = base ** exponent

# Print the final equation and the result, showing each number.
print(f"The graph is K_1,{n}.")
print(f"The global labeling number is the {n}th term of the sequence x_k = 2^(k-1).")
print(f"The final equation is: x_{n} = {base}^({n}-1) = {base}^{exponent}")
print(f"The calculated global labeling number is: {global_labeling_number}")