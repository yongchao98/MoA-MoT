import math

# Step 1: Define the given parameters from the problem description.

# The frequency of the homozygous recessive genotype (0/0), which cannot digest milk.
# The problem states this is half the population.
freq_no_milk_provider = 0.5

# The frequency of individuals who can digest milk is the rest of the population.
freq_milk_provider = 1 - freq_no_milk_provider

# Height for individuals who had milk during development, converted to inches.
# 4 feet 6 inches = (4 * 12) + 6
height_with_milk_in = (4 * 12) + 6

# Height for individuals without milk during development, converted to inches.
# 3 feet 6 inches = (3 * 12) + 6
height_without_milk_in = (3 * 12) + 6

# Step 2: Calculate the components of the weighted average.
# The children's height depends on whether their father can provide milk.
# The proportion of children who get milk is equal to the frequency of fathers who can provide it.
prop_with_milk = freq_milk_provider
prop_without_milk = freq_no_milk_provider

# Calculate the contribution of each group to the average height.
height_contribution_with_milk = prop_with_milk * height_with_milk_in
height_contribution_without_milk = prop_without_milk * height_without_milk_in

# Step 3: Calculate the final average height.
average_height = height_contribution_with_milk + height_contribution_without_milk

# Step 4: Print the full calculation and the result formatted to four significant figures.
print("Calculating the average adult height of the population.")
print("-" * 50)
print("The proportion of the population that grows up with milk is equal to the proportion of fathers who can digest lactose.")
print(f"Proportion with milk = 1 - Frequency(0/0) = 1 - {prop_without_milk} = {prop_with_milk}")
print("\nThe proportion of the population that grows up without milk is equal to the proportion of fathers who cannot digest lactose.")
print(f"Proportion without milk = Frequency(0/0) = {prop_without_milk}")
print("\nHeights in inches:")
print(f"Height with milk (4'6\") = {height_with_milk_in} inches")
print(f"Height without milk (3'6\") = {height_without_milk_in} inches")
print("\nThe average height is the weighted average:")
print("Average Height = (Proportion with milk * Height with milk) + (Proportion without milk * Height without milk)")
# Print the equation with the actual numbers
print(f"Average Height = ({prop_with_milk} * {height_with_milk_in}) + ({prop_without_milk} * {height_without_milk_in})")
print(f"Average Height = {height_contribution_with_milk} + {height_contribution_without_milk}")
# The result is 48, which has two significant figures. To show four, we format it as 48.00.
print(f"Final Average Height = {average_height:.2f} inches")
print("-" * 50)

# The final answer in the required format. The calculated value is 48.
# To round to four significant figures, we present it as 48.00.
final_answer_val = f"{average_height:.2f}"
# The final output tag is <<<48.00>>>