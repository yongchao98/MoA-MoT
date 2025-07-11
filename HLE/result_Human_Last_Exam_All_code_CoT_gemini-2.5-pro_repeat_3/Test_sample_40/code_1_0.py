import math

# Step 1: Define the proportions of the population with and without milk.
# The proportion of the population without milk is determined by the frequency
# of fathers with the 0/0 genotype, which is given as 0.5.
prop_without_milk = 0.5
prop_with_milk = 1 - prop_without_milk

# Step 2: Define the heights in inches for each group.
# Height without milk is 3 feet 6 inches.
height_without_milk = 3 * 12 + 6
# Height with milk is 1 foot taller, so 4 feet 6 inches.
height_with_milk = height_without_milk + 12

# Step 3: Calculate the average height of the population.
# This is the weighted average of the heights of the two groups.
average_height = (prop_without_milk * height_without_milk) + (prop_with_milk * height_with_milk)

# Step 4: Print the final equation with all the numbers.
print("The average height is calculated using the following equation:")
print(f"({prop_without_milk} * {height_without_milk}) + ({prop_with_milk} * {height_with_milk}) = {average_height}")

# Format the final answer to four significant figures. 48.00 has four.
final_answer_formatted = f"{average_height:.2f}"

print(f"\nThe average adult height of the population is {final_answer_formatted} inches.")

# The final numerical answer for the submission format.
# Note: This print is for internal use to generate the final answer tag.
# print(f"\n<<<...>>>\n{final_answer_formatted}")