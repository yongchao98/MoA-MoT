import math

# The optimal success probability for Agent C, p_C, is derived to be 0.5.
# This value is the supremum of C's winning probability function.
# C can achieve a winning probability arbitrarily close to 0.5 by choosing
# a value just under 0.5.
p_C = 0.5

# The problem asks for the floor of 100 * p_C.
result = math.floor(100 * p_C)

# Print the final equation and the result
print(f"p_C = {p_C}")
print(f"floor(100 * p_C) = floor(100 * {p_C}) = floor({100 * p_C}) = {result}")
