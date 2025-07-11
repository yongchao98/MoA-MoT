# Step 1: Define the proportions of the population for each height group.
# The proportion of fathers with the 0/0 genotype determines the proportion of children without milk.
prop_no_milk = 0.5
# The remaining proportion of fathers can provide milk.
prop_with_milk = 1 - prop_no_milk

# Step 2: Define the heights in inches.
# Height for individuals without milk during development (3 feet 6 inches).
height_no_milk_inches = 3 * 12 + 6
# Height for individuals with milk during development (4 feet 6 inches).
height_with_milk_inches = 4 * 12 + 6

# Step 3: Calculate the weighted average height of the population.
# Average height = (proportion without milk * height without milk) + (proportion with milk * height with milk)
average_height = (prop_no_milk * height_no_milk_inches) + (prop_with_milk * height_with_milk_inches)

# Step 4: Print the final equation and the result.
# The problem asks to output each number in the final equation.
print("The average adult height of the population is calculated by taking the weighted average of the two possible heights.")
print(f"Calculation: ({prop_no_milk} * {height_no_milk_inches}) + ({prop_with_milk} * {height_with_milk_inches}) = {average_height:.4g}")
print(f"The average adult height is {average_height:.4g} inches.")
# Final answer formatted to four significant figures as 48.00.
print(f"The average adult height is {average_height:.2f} inches.")

# Return the answer in the specified format. The final value is 48, which needs to be represented with four significant figures.
final_answer_str = f"{average_height:.2f}"
# The wrapper format is not required here, but we will print the final numerical answer.
# For example: <<<48.00>>>
# As per instructions, printing the final result at the end is what matters.
