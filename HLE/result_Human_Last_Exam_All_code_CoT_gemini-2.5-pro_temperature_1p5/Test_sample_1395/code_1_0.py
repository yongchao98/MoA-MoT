# Step 1: Determine the chromatic number of the Asian country graph before 1991.
# A map graph is planar, so its chromatic number is at most 4.
# A 4-coloring is required if there is a country surrounded by an odd-number of mutually bordering countries.
# Before 1991, Laos was bordered by five countries (China, Myanmar, Thailand, Cambodia, Vietnam) that form a 5-cycle.
# This forms a W6 wheel graph, which requires 4 colors.
chromatic_number_before = 4
print(f"The chromatic number of the Asian graph before the Soviet dissolution was {chromatic_number_before}.")

# Step 2: Determine the chromatic number of the Asian country graph after 1991.
# After the dissolution, the structure around Laos still exists, keeping the chromatic number at 4.
# The creation of new Central Asian republics also creates a new 4-color-forcing structure, but does not raise the overall number.
chromatic_number_after = 4
print(f"The chromatic number of the Asian graph after the Soviet dissolution is {chromatic_number_after}.")

# Step 3: Calculate delta_soviet, the change in the chromatic number.
delta_soviet = chromatic_number_after - chromatic_number_before
print(f"The incremental change, delta_soviet, is {chromatic_number_after} - {chromatic_number_before} = {delta_soviet}.")

# Step 4: Determine beta, representing the change in planarity.
# A map of contiguous countries is a planar graph. Subdividing a country preserves planarity.
# The problem states beta = 1 if planarity does not change.
beta = 1
print(f"The dissolution did not change the graph's planarity, so beta = {beta}.")

# Step 5: Calculate and print the final result.
final_answer = beta * delta_soviet
print(f"The final result is beta * delta_soviet, which is calculated as:")
print(f"{beta} * {delta_soviet} = {final_answer}")