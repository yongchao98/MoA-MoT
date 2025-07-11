# Step 1: Define the proportions and heights based on the problem description.

# The problem states half the population has the 0/0 genotype.
# A family has milk only if the father does NOT have the 0/0 genotype.
# Therefore, the proportion of the population that grows up without milk is 0.5.
prop_without_milk = 0.5

# The remaining proportion of the population grows up with milk.
prop_with_milk = 1 - prop_without_milk

# Step 2: Convert heights from feet and inches to total inches.
# Height with milk is 4 feet 6 inches.
height_with_milk_inches = 4 * 12 + 6

# Height without milk is 3 feet 6 inches.
height_without_milk_inches = 3 * 12 + 6

# Step 3: Calculate the average height of the population.
# This is the weighted average of the two heights.
avg_height = (prop_with_milk * height_with_milk_inches) + (prop_without_milk * height_without_milk_inches)

# Step 4: Print the final equation and the result.
# The problem requests that each number in the final equation be output.
print(f"The average height is calculated as: ({prop_with_milk} * {height_with_milk_inches}) + ({prop_without_milk} * {height_without_milk_inches})")

# The result needs to be formatted to four significant figures.
# The calculated value is 48. To represent this with four significant figures, we format it as 48.00.
print(f"The average adult height of the population is {avg_height:.2f} inches.")

# Print the final answer in the required format
print(f"<<<{avg_height:.2f}>>>")