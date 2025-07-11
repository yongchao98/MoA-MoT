# 1. Determine the chromatic number of the Asian country graph before the Soviet Union's dissolution.
# A map of countries is a planar graph. By the Four-Color Theorem, its chromatic number is at most 4.
# The complex pre-dissolution map of Asia required 4 colors.
chromatic_number_before = 4
print(f"Chromatic number of the Asian graph before dissolution = {chromatic_number_before}")

# 2. Determine the chromatic number after the dissolution.
# The dissolution increased the number of countries and borders, making the graph more complex but still planar.
# The graph still requires 4 colors, as no simplification occurred.
chromatic_number_after = 4
print(f"Chromatic number of the Asian graph after dissolution = {chromatic_number_after}")

# 3. Calculate the incremental change, delta_soviet.
delta_soviet = chromatic_number_after - chromatic_number_before
print(f"The change, delta_soviet, is {chromatic_number_after} - {chromatic_number_before} = {delta_soviet}")

# 4. Determine the change in planarity for beta.
# A political map is inherently planar. Subdividing a country preserves the graph's planarity.
# Since the planarity did not change, beta is 1.
beta = 1
print(f"Planarity did not change, so beta = {beta}")

# 5. Calculate the final result.
final_answer = beta * delta_soviet
print(f"The final result is {beta} * {delta_soviet} = {final_answer}")