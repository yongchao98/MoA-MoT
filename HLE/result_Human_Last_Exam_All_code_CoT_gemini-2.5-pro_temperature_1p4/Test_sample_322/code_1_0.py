# The Three Utilities Problem: A Mathematical Check

# 1. Define the graph K(3,3)
# There are 3 houses and 3 utilities.
num_houses = 3
num_utilities = 3

# V is the total number of vertices (nodes).
V = num_houses + num_utilities

# E is the total number of edges (utility lines).
# Each house must connect to each of the 3 utilities.
E = num_houses * num_utilities

print(f"Analyzing the Three Utilities Problem (Graph K3,3):")
print(f"Number of vertices (V) = {V}")
print(f"Number of edges (E) = {E}")
print("-" * 30)

# 2. State the condition for planarity in a bipartite graph.
# For any simple, connected, bipartite planar graph, the inequality E <= 2V - 4 must hold true.
# This is a consequence of Euler's formula for planar graphs (V - E + F = 2) combined
# with the fact that bipartite graphs have no cycles of odd length (i.e., the smallest
# possible face is bounded by 4 edges).

print("Checking the planarity condition for bipartite graphs: E <= 2*V - 4")
print("-" * 30)

# 3. Substitute the values of V and E and check the inequality.
inequality_right_side = 2 * V - 4
is_planar = E <= inequality_right_side

# Print the final equation with the numbers substituted in.
print(f"Substituting the values for K3,3:")
print(f"{E} <= 2 * {V} - 4")
print(f"{E} <= {2 * V} - 4")
print(f"{E} <= {inequality_right_side}")
print("-" * 30)


# 4. State the conclusion.
print(f"Is the condition met? {is_planar}")

if not is_planar:
    print("\nThe inequality is false. This provides a mathematical proof that the graph K3,3 is non-planar.")
    print("Therefore, it is impossible to draw it on a 2D plane without at least one crossing.")
    print("This means the puzzle, under the strict rules provided, has no solution.")
else:
    print("\nThe inequality is true, suggesting the graph could be planar.")
    print("However, the check is a necessary but not sufficient condition. Kuratowski's theorem is the definitive proof.")

<<<E>>>