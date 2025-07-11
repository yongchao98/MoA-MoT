# Based on the topological proof, the problem is to determine a specific integer value.
# The proof proceeds in two parts:
# 1. Establish a lower bound for the cardinality.
# 2. Construct an example that meets this lower bound.

# Part 1: Lower Bound
# For any decomposable continuum X = A U B, we can derive at least two
# distinct regular proper subcontinua from the components of X \ A and X \ B.
# This proves the cardinality must be at least 2.
lower_bound = 2

# Part 2: Achievable Example
# An example can be constructed by joining two pseudo-arcs at a single point.
# This continuum is decomposable, and its only two regular proper subcontinua
# are the pseudo-arcs themselves.
# This shows a cardinality of 2 is possible.
example_cardinality = 2

# The smallest possible cardinality is the minimum of all possible cardinalities,
# which is established by the lower bound that is shown to be achievable.
smallest_possible_cardinality = min(n for n in [example_cardinality] if n >= lower_bound)

# The final answer is a single number.
print(smallest_possible_cardinality)