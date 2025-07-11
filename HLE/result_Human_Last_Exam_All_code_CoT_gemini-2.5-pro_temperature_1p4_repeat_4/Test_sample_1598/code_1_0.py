# The graph is K_1,100, which has 100 edges.
# The problem of finding the global labeling number reduces to finding the
# 100th term of a sequence of labels {w_1, w_2, ..., w_100} such that
# for any k, w_k is not a subset sum of the other labels.
# The optimal sequence that minimizes the maximum label is the sequence of
# powers of two: w_k = 2**(k-1).
# We need to find the 100th term of this sequence.

# The number of edges in the graph determines the size of the label set.
num_edges = 100

# The equation for the k-th label is w_k = 2^(k-1).
# For the 100th label, the equation is 2^(100-1) = 2^99.
base = 2
exponent = num_edges - 1

# Calculate the result
result = base ** exponent

# The global labeling number is the largest label required, which is the final result.
# The final equation is 2^99. Let's output all numbers in this equation.
print(f"The global labeling number of K_1,100 is given by the equation: {base}^{exponent}")
print("This evaluates to:")
print(result)