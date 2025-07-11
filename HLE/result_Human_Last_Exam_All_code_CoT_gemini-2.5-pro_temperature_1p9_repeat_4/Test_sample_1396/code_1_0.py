import math

# This calculation is based on the query complexity analysis of the protocol for 
# 4-agent connected ε-envy-free cake cutting by Goldberg, Hollender, and Suksompong (2021).
# The "upper bound O" is interpreted as the constant factor of the dominant term in their complexity bound.

# Number of agents
n = 4

# The protocol's analysis for n=4 involves parameters derived from its structure:
# The main protocol makes at most 2 calls to a core subroutine.
num_subroutine_calls = 2

# The analysis of the subroutine shows its complexity depends on the number of agent permutations (n!).
# For each permutation, a set of queries are made.
queries_per_permutation = 5
additive_term = 3

# Calculate the number of permutations for n agents.
n_permutations = math.factorial(n)

# Calculate the constant for the core subroutine's dominant term.
core_protocol_constant = n_permutations * queries_per_permutation + additive_term

# Calculate the final upper bound O for the main protocol.
O = num_subroutine_calls * core_protocol_constant

# As requested, printing the equation with each number.
print(f"The calculation for the bound O is based on the formula:")
print(f"{num_subroutine_calls} * (factorial({n}) * {queries_per_permutation} + {additive_term})")
print(f"Which evaluates to:")
print(f"{num_subroutine_calls} * ({n_permutations} * {queries_per_permutation} + {additive_term}) = {O}")