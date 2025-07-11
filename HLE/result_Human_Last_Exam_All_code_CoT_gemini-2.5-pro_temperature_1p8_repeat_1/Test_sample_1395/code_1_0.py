import sys

# This script calculates a value based on changes to the graph of Asian countries
# following the dissolution of the Soviet Union.

# Step 1: Determine the chromatic number of the Asian graph AFTER the dissolution.
# Post-1991, new Central Asian states were formed. Let's analyze the subgraph around Uzbekistan.
# Uzbekistan is bordered by: Kazakhstan, Kyrgyzstan, Tajikistan, Afghanistan, and Turkmenistan.
# These five neighbors form a 5-cycle border arrangement:
# Kazakhstan -> Kyrgyzstan -> Tajikistan -> Afghanistan -> Turkmenistan -> Kazakhstan.
# A 5-cycle (an odd cycle) requires 3 colors to be colored correctly.
# Uzbekistan, as a central node, borders all five of these countries.
# The colors used for the 5-cycle are, for example, C1, C2, C1, C2, C3.
# Uzbekistan is adjacent to nodes with colors C1, C2, and C3. Therefore, it needs a fourth color, C4.
# This structure (a W6 wheel graph) forces the graph to be 4-chromatic.
chromatic_number_after = 4
print(f"Chromatic number of the Asian graph after Soviet dissolution: {chromatic_number_after}")

# Step 2: Determine the chromatic number of the Asian graph BEFORE the dissolution.
# Before 1991, the region of Uzbekistan and its neighbors was a single entity: the USSR.
# The W6 wheel graph structure described above did not exist.
# While large political maps can be complex, the most direct consequence of the dissolution
# was the creation of this specific 4-chromatic structure. We can infer that before
# this change, the graph was 3-colorable.
chromatic_number_before = 3
print(f"Chromatic number of the Asian graph before Soviet dissolution: {chromatic_number_before}")

# Step 3: Calculate delta_soviet, the incremental change.
delta_soviet = chromatic_number_after - chromatic_number_before
print(f"The incremental change, delta_soviet, is {chromatic_number_after} - {chromatic_number_before} = {delta_soviet}")

# Step 4: Determine beta based on the change in planarity.
# A graph of bordering sovereign nations on a map is, by definition, a planar graph.
# The dissolution of the USSR created new nodes and edges, but the resulting graph is still
# a map of countries and remains planar.
# The problem states if planarity does not change, beta = 1.
beta = 1
print(f"The planarity of the graph did not change, so beta = {beta}")

# Step 5: Calculate the final result.
final_answer = beta * delta_soviet
print(f"The final answer is beta * delta_soviet = {beta} * {delta_soviet} = {final_answer}")

# Output the final answer in the required format.
# The triple angle brackets are for the system to parse the answer.
sys.stdout.flush()
print(f"<<<{final_answer}>>>")