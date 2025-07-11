import sys

# This script calculates a value based on changes in the graph theory properties of the Asian political map
# due to the dissolution of the Soviet Union.

# Step 1: Determine the chromatic number of the Asian map BEFORE the dissolution.
# Before 1991, the USSR was a single large nation. The political map of Asia was complex, but it is
# generally understood in graph theory that this map was 3-colorable. The graph contains 3-cliques
# (e.g., USSR, China, Afghanistan), so its chromatic number is at least 3. A 3-coloring is possible.
chi_before = 3
print(f"Chromatic number of the Asian country graph before dissolution (chi_before): {chi_before}")

# Step 2: Determine the chromatic number of the Asian map AFTER the dissolution.
# The breakup of the USSR created several new countries in Central Asia, significantly increasing the
# density of borders in that region. This new, more complex map is no longer 3-colorable and is
# known to require 4 colors, which is the maximum for any planar graph by the Four Color Theorem.
chi_after = 4
print(f"Chromatic number of the Asian country graph after dissolution (chi_after): {chi_after}")

# Step 3: Calculate the incremental change, delta_soviet.
# The change is the difference between the new chromatic number and the old one.
delta_soviet = chi_after - chi_before
print(f"The incremental change, delta_soviet = chi_after - chi_before = {chi_after} - {chi_before} = {delta_soviet}")

# Step 4: Determine if the planarity of the graph changed.
# A political map of contiguous nations on a globe is inherently a planar graph.
# The dissolution of one country into several new ones does not change this fundamental property.
# The graph remained planar. According to the problem, if planarity does not change, beta = 1.
beta = 1
print(f"The graph's planarity did not change, so beta = {beta}")

# Step 5: Calculate and print the final answer.
final_answer = beta * delta_soviet
print(f"Final Answer = beta * delta_soviet = {beta} * {delta_soviet} = {final_answer}")

# Final answer in the required format for the system.
sys.stdout.write(f"\n<<<{final_answer}>>>\n")
