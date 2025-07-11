# Step 1: Determine the chromatic number of the Asian map graph before the Soviet dissolution.
# A map of countries is a planar graph. The Four Color Theorem guarantees its chromatic number is at most 4.
# The complexity of the Asian map requires 4 colors.
chromatic_number_before = 4
print(f"The chromatic number of the Asian map graph before the dissolution was {chromatic_number_before}.")

# Step 2: Determine the chromatic number of the Asian map graph after the Soviet dissolution.
# After the dissolution, new countries were formed, but the map remains a planar graph.
# The Four Color Theorem still applies, and the chromatic number is still 4.
chromatic_number_after = 4
print(f"The chromatic number of the Asian map graph after the dissolution is {chromatic_number_after}.")

# Step 3: Calculate delta_soviet, the incremental change in the chromatic number.
delta_soviet = chromatic_number_after - chromatic_number_before
print(f"The change in chromatic number, delta_soviet, is {chromatic_number_after} - {chromatic_number_before} = {delta_soviet}.")

# Step 4: Determine if the planarity of the graph changed and assign a value to beta.
# A country map is inherently a planar graph. Adding new countries and borders does not change this property.
# The problem states that if planarity does not change, beta = 1.
beta = 1
print(f"The planarity of the graph did not change, so beta = {beta}.")

# Step 5: Calculate and print the final result.
final_answer = beta * delta_soviet
print(f"The final result is beta * delta_soviet = {beta} * {delta_soviet} = {final_answer}.")