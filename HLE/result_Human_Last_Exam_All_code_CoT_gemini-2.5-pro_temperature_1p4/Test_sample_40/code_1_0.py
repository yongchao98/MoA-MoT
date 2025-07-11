# Step 1: Define the proportions of the population with and without milk.
# The problem states that half of the population has the 0/0 genotype.
# Since milk availability depends on the father's genotype and mating is random,
# the proportion of families without milk is the same as the frequency of the 0/0 genotype.
prop_without_milk = 0.5
prop_with_milk = 1 - prop_without_milk

# Step 2: Define the heights in inches.
# Height without milk: 3 feet 6 inches = (3 * 12) + 6 inches
height_without_milk_in = 42
# Height with milk: 1 foot taller than without milk = 4 feet 6 inches = (4 * 12) + 6 inches
height_with_milk_in = 54

# Step 3: Calculate the weighted average height of the population.
# Average Height = (proportion with milk * height with milk) + (proportion without milk * height without milk)
average_height = (prop_with_milk * height_with_milk_in) + (prop_without_milk * height_without_milk_in)

# Step 4: Print the final equation and the result, formatted to four significant figures.
# The result is 48. To express it with four significant figures, we write 48.00.
# The format specifier ':.2f' will achieve this for the number 48.
print("The average adult height is the weighted average of the heights of individuals who grew up with milk and those who grew up without.")
print(f"Final Equation: ({prop_with_milk} * {height_with_milk_in}) + ({prop_without_milk} * {height_without_milk_in}) = {average_height:.2f}")
print(f"The average adult height of the population is {average_height:.2f} inches.")
