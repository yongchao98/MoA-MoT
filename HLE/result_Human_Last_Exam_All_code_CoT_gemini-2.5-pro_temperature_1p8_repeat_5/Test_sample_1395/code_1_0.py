# The chromatic number of a graph is the minimum number of colors needed to color its vertices
# so that no two adjacent vertices share the same color. A graph of nations with shared borders is a planar graph.

# By the Four-Color Theorem, the chromatic number of any planar graph is at most 4.
# The complex web of borders in Asia before 1991 makes it a standard case requiring 4 colors.
chi_before = 4
print(f"The chromatic number of the Asian graph before the Soviet dissolution, chi_before: {chi_before}")

# After the dissolution, the USSR node was replaced by multiple new state nodes.
# This made the graph more complex, but it remained planar. The pre-existing complexity that
# required 4 colors (e.g., in other parts of Asia) was not simplified. Thus, 4 colors are still necessary.
chi_after = 4
print(f"The chromatic number of the Asian graph after the Soviet dissolution, chi_after: {chi_after}")

# delta_soviet is the incremental change in the chromatic number.
delta_soviet = chi_after - chi_before
print(f"The incremental change, delta_soviet = chi_after - chi_before = {chi_after} - {chi_before} = {delta_soviet}")

# A graph of nations on a map is planar. The dissolution of one nation into several
# new nations is an operation that preserves the planarity of the graph.
# Therefore, the planarity did not change, so beta is 1.
beta = 1
print(f"The planarity of the graph did not change, so beta: {beta}")

# The final result is the product of beta and delta_soviet.
final_answer = beta * delta_soviet
print(f"The final result is beta * delta_soviet = {beta} * {delta_soviet} = {final_answer}")

print(f"\n<<<0>>>")