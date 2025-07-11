# Step 1: Define the chromatic number of the Asian country graph before the Soviet dissolution.
# Before 1991, the map of Asia was 3-colorable because no four countries were all mutually adjacent.
chi_before = 3
print(f"The chromatic number of the Asian country graph before the dissolution was: {chi_before}")

# Step 2: Define the chromatic number after the Soviet dissolution.
# The dissolution created Tajikistan, which borders China, Afghanistan, and Pakistan. 
# These four countries form a 4-clique, requiring a minimum of 4 colors.
chi_after = 4
print(f"The chromatic number of the Asian country graph after the dissolution is: {chi_after}")

# Step 3: Calculate the incremental change, delta_soviet.
delta_soviet = chi_after - chi_before
print(f"The incremental change, delta_soviet, is the difference: {chi_after} - {chi_before} = {delta_soviet}")

# Step 4: Determine the value of beta based on planarity.
# A map graph is always planar. The political change of subdividing a country does not alter this property.
# Since the planarity of the graph did not change, beta is 1.
beta = 1
print(f"The planarity of the graph did not change, therefore beta is: {beta}")

# Step 5: Calculate and print the final result.
final_answer = beta * delta_soviet
print(f"The final answer is beta * delta_soviet = {beta} * {delta_soviet} = {final_answer}")