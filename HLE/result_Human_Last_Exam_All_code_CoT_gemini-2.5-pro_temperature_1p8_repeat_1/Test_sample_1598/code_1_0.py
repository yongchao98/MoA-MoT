# The graph is K_1,100, which has n=100 edges from a central vertex.
num_edges = 100

# As deduced from the problem's definition, the global labeling condition for K_1,n
# requires a set of n labels where no label can be represented as a sum of a
# subset of the other labels.
#
# To minimize the maximum label (the global labeling number), we must choose the
# smallest possible integers that satisfy this property. This leads to a sequence
# where each new label is greater than the sum of all preceding labels.
#
# The optimal sequence starting with the smallest possible integer (1) is:
# w_1 = 1
# w_2 = 2
# w_3 = 4
# ...
# w_i = 2**(i-1)
#
# The global labeling number is the maximum label needed, which is the 100th
# label in this sequence, w_100.

# The base of the power.
base = 2

# The exponent for the 100th term (since the sequence starts with 2**0).
exponent = num_edges - 1

# Calculate the result.
result = base ** exponent

print(f"The global labeling number of K_1,{num_edges} is determined by the maximum value in an optimally chosen set of labels.")
print("This corresponds to the 100th term in the sequence of powers of two.")
print(f"The final equation to solve is: {base} ** ({num_edges} - 1)")
print(f"Result: {result}")