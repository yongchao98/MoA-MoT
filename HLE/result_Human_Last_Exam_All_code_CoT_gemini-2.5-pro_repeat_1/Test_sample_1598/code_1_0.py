import math

# For the graph K_{1,100}, the problem of finding the global labeling number
# reduces to finding a set of 100 positive integers {w_1, ..., w_100} such that
# no element is a sum of any subset of the other elements, and the maximum
# element in the set is minimized.

# This is achieved by constructing the set greedily. Let the labels be sorted:
# w_1 < w_2 < ... < w_100.
# The condition means that for any i, w_i cannot be a sum of a subset of {w_1, ..., w_{i-1}}.

# The sequence that minimizes the labels is constructed as follows:
# w_1 = 1
# w_2 = 2 (since it can't be 1)
# w_3 = 4 (since it can't be 1, 2, or 1+2=3)
# In general, the i-th label in this optimal sequence is w_i = 2^(i-1).

# We need to find the 100th label for the graph K_{1,100}.
n = 100
base = 2
power = n - 1

# Calculate the result, which is the global labeling number.
# Python's integers handle arbitrary size, so this calculation is exact.
result = base ** power

# The problem asks us to determine the global labeling number.
# The following lines print the logic and the final calculated number.
print(f"The global labeling number for K_1,100 is determined by the 100th element in the sequence w_i = 2^(i-1).")
print(f"This requires calculating the value of the final and largest label, w_100.")
print(f"The final equation is: {base} ** {power}")
print(f"Result: {result}")