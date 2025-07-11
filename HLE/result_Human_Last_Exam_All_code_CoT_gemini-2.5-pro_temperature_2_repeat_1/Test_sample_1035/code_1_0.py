import math

# The number of vertices in the complete graph K_n
n = 35

# According to a well-known theorem by Graham and Pollak, the minimum number
# of bipartite graphs needed to cover all the edges of a complete graph K_n
# is given by the ceiling of the base-2 logarithm of n.
# This is also known as the bipartite dimension of K_n.

# The formula is k = ceil(log2(n))

# We calculate the base-2 logarithm of n.
log_val = math.log2(n)

# Then we take the ceiling of this value to get the minimum number of graphs.
result = math.ceil(log_val)

print(f"The complete graph is K_n where n = {n}.")
print(f"The problem is to find the minimum number of bipartite graphs needed to cover the edges of K_{n}.")
print(f"This is calculated by the formula: ceil(log2(n)).")
print(f"Step 1: Calculate log2(n) for n = {n}")
print(f"log2({n}) = {log_val}")
print(f"Step 2: Calculate the ceiling of the result.")
print(f"ceil({log_val}) = {result}")
print(f"Final Answer: The minimum number of bipartite graphs is {result}.")