# 1. Define the known variables based on the problem description.

# The proportion of the population with genotype 0/0. These are the individuals
# who cannot digest milk. Since milk availability depends on the father's genotype
# and mating is random, this is also the proportion of children whose father is 0/0.
proportion_no_milk_father = 0.5

# The proportion of children whose father can provide milk (genotype 1/0 or 1/1).
# This is the rest of the population.
proportion_with_milk_father = 1.0 - proportion_no_milk_father

# Convert the heights from feet and inches to total inches.
# Tall height (with milk): 4 feet * 12 inches/foot + 6 inches
height_with_milk = 4 * 12 + 6

# Short height (without milk): 3 feet * 12 inches/foot + 6 inches
height_no_milk = 3 * 12 + 6

# 2. Calculate the average height of the population.
# This is a weighted average based on the proportions of tall and short individuals.
average_height = (proportion_no_milk_father * height_no_milk) + (proportion_with_milk_father * height_with_milk)

# 3. Print the equation and the final result.
# The problem asks to output each number in the final equation.

print("The average height is a weighted average of the two possible heights in the population.")
print("Average Height = (Proportion without milk * Height without milk) + (Proportion with milk * Height with milk)")
print(f"Average Height = ({proportion_no_milk_father} * {height_no_milk}) + ({proportion_with_milk_father} * {height_with_milk})")
print(f"Average Height = {proportion_no_milk_father * height_no_milk} + {proportion_with_milk_father * height_with_milk}")
# The final result must be formatted to four significant figures.
# The exact answer is 48. To show four significant figures, we format it as 48.00.
print(f"Average Height = {average_height:.2f} inches")

print("\n<<<48.00>>>")