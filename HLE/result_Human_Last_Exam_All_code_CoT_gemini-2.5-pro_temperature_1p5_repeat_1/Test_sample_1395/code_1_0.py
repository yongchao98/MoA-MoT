# Step 1: Calculate the change in the chromatic number (delta_soviet)

# The chromatic number of the complex Asian country map before the dissolution of the Soviet Union.
# Based on the Four-Color Theorem and the map's complexity, it required 4 colors.
chromatic_number_before = 4
print(f"Chromatic number of the Asian subgraph before the dissolution: {chromatic_number_before}")

# The chromatic number after the dissolution. The map is still planar (max 4 colors)
# and contains the previous complexity, so it also requires 4 colors.
chromatic_number_after = 4
print(f"Chromatic number of the Asian subgraph after the dissolution: {chromatic_number_after}")

# Calculating the change.
delta_soviet = chromatic_number_after - chromatic_number_before
print(f"The incremental change, delta_soviet, is {chromatic_number_after} - {chromatic_number_before} = {delta_soviet}")

# Step 2: Determine beta based on the change in planarity.

# A map of countries is planar by definition.
is_planar_before = True
is_planar_after = True

# Check if the planarity status changed.
planarity_changed = (is_planar_before != is_planar_after)

# If planarity did not change, beta is 1. If it changed, beta is 2.
if planarity_changed:
    beta = 2
else:
    beta = 1
print(f"\nThe planarity of the graph did not change, therefore beta = {beta}")

# Step 3: Calculate the final answer.
final_answer = beta * delta_soviet
print(f"\nThe final result is beta * delta_soviet = {beta} * {delta_soviet} = {final_answer}")