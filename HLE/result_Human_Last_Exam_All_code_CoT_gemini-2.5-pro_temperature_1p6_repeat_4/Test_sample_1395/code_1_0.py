# Step 1: Determine the chromatic number of the Asian subgraph before the Soviet dissolution.
# A graph of countries on a map is planar, so its chromatic number is at most 4 by the Four-Color Theorem.
# To determine the exact number, consider the subgraph formed by the USSR, China, Afghanistan, and Pakistan.
# A valid 4-coloring can be forced as follows:
# 1. Assign China 'Color 1'.
# 2. The USSR borders China, so it must be a different color. Assign USSR 'Color 2'.
# 3. Pakistan borders both China ('Color 1') and the USSR ('Color 2'). It requires 'Color 3'.
# 4. Afghanistan borders China ('Color 1'), the USSR ('Color 2'), and Pakistan ('Color 3').
# Therefore, Afghanistan needs a fourth distinct color, 'Color 4'.
# This proves that 4 colors were necessary.
chi_before = 4
print(f"The chromatic number of the Asian graph before the dissolution was: {chi_before}")

# Step 2: Determine the chromatic number after the Soviet dissolution.
# The graph remains planar, so the chromatic number is still at most 4.
# The dissolution created new countries in Central Asia. Consider Uzbekistan and its five neighbors:
# Kazakhstan, Kyrgyzstan, Tajikistan, Afghanistan, and Turkmenistan.
# These five neighbors form a 5-cycle on the map (an odd-length cycle), which requires 3 colors.
# Uzbekistan borders all five of these countries in the cycle.
# To color Uzbekistan, a fourth color is required, as it's adjacent to nodes of all three colors used for the cycle.
# This structure (a central node connected to an odd cycle) proves the chromatic number is 4.
chi_after = 4
print(f"The chromatic number of the Asian graph after the dissolution is: {chi_after}")

# Step 3: Calculate the change in chromatic number, delta_soviet.
# delta_soviet is the difference between the new chromatic number and the old one.
delta_soviet = chi_after - chi_before
print(f"The change in the chromatic number, delta_soviet = {chi_after} - {chi_before} = {delta_soviet}")

# Step 4: Determine the change in planarity.
# A graph based on a geographical map is planar. Subdividing a region (the USSR) into smaller regions (new countries)
# does not introduce any line crossings that would violate planarity.
# The planarity of the graph did not change.
# According to the problem, beta is 1 if planarity does not change.
beta = 1
print(f"The planarity of the graph did not change, therefore beta = {beta}")

# Step 5: Calculate the final result.
final_answer = beta * delta_soviet
print(f"The final calculation is beta * delta_soviet = {beta} * {delta_soviet} = {final_answer}")
print(f"<<<{final_answer}>>>")