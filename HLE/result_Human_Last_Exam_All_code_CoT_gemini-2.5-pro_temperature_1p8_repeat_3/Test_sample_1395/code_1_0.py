# Step 1: Determine the chromatic number of the Asian graph before the Soviet dissolution.
# A political map is a planar graph. By the Four-Color Theorem, the chromatic number is at most 4.
# Due to the high complexity of borders in Asia, the graph required 4 colors.
chromatic_number_before = 4
print(f"Chromatic number of the Asian graph before 1991 (χ_before): {chromatic_number_before}")

# Step 2: Determine the chromatic number after the Soviet dissolution.
# The dissolution added new states and borders but did not fundamentally change the graph's nature.
# It remained a complex planar graph requiring 4 colors.
chromatic_number_after = 4
print(f"Chromatic number of the Asian graph after 1991 (χ_after): {chromatic_number_after}")

# Step 3: Calculate the incremental change, δ_soviet.
delta_soviet = chromatic_number_after - chromatic_number_before
print(f"The incremental change δ_soviet = {chromatic_number_after} - {chromatic_number_before} = {delta_soviet}")

# Step 4: Determine the change in planarity and the value of β.
# The map of countries was and remains a planar graph. Subdividing a region does not change planarity.
# The problem states β = 1 if planarity does not change.
planarity_changed = False
if not planarity_changed:
    beta = 1
else:
    beta = 2
print(f"The planarity of the graph did not change, so β = {beta}")

# Step 5: Calculate the final result.
final_result = beta * delta_soviet
print(f"Final Answer (β * δ_soviet) = {beta} * {delta_soviet} = {final_result}")
