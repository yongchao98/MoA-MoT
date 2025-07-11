import math

# Step 1: Calculate the external arrangement factor.
# This involves arranging the CSM-Block, 2 Ethicists, and 4 other Classicists.
# We arrange the CSM-Block and 2 Ethicists in a circle: (3-1)! ways.
# The 4 other Classicists must sit in the single slot between the Ethicists,
# and can be arranged in 4! ways.
external_arrangements = math.factorial(3 - 1) * math.factorial(4)

# Step 2: Calculate the internal arrangement factor for the CSM-Block.
# From our casework, this is 152 multiplied by the arrangements of the 9 non-end, non-rower scientists (9!).
# The full derivation is: (ways_S_Female_M_Male * 1_spot) + (ways_S_Male_M_Female * 1_spot) + (ways_S_Female_M_Female * 2_spots)
# This simplifies to (24 + 32 + 96) * 9! = 152 * 9!
internal_arrangement_factor = 152

# Step 3: Calculate the factorial term.
factorial_term = math.factorial(9)

# Step 4: Calculate the total number of arrangements.
total_ways = external_arrangements * internal_arrangement_factor * factorial_term

# Final Output: Print the components of the final equation as requested.
print(f"The final calculation is the product of external arrangements, the internal arrangement factor, and the factorial term.")
print(f"{external_arrangements} * {internal_arrangement_factor} * {factorial_term} = {total_ways}")
