# The problem is to find the minimum length n of an observation sequence
# such that a 3-state memory can outperform a 2-state memory, which in
# turn is no better than a 1-state (memoryless) agent.

# This corresponds to finding the length of the shortest binary strings
# that can be distinguished by a 3-state automaton (k=3) but not by any
# 2-state automaton (k-1=2).

# The size of the automaton that can distinguish the sequences is k.
k = 3

# The formula for the minimum length n (often denoted g(k) in automata theory)
# for an alphabet of size 2 (our observations are binary) is:
# n = 2^(k-1) + k - 2

# We substitute k=3 into the formula.
term1 = 2**(k - 1)
term2 = k - 2
n = term1 + term2

print(f"The minimum hallway length 'n' is calculated using the formula for string distinguishability by automata.")
print(f"Let k be the number of memory states that allows for a higher return, so k = {k}.")
print(f"The formula is n = 2^(k-1) + k - 2.")
print(f"Substituting k = {k}, the equation becomes:")
# In the final response, remember to output each number in the final equation.
print(f"n = 2^({k} - 1) + {k} - 2")
print(f"n = 2^{k-1} + {term2}")
print(f"n = {term1} + {term2}")
print(f"n = {n}")

print("\nTherefore, the minimum length of the hallway is 5.")