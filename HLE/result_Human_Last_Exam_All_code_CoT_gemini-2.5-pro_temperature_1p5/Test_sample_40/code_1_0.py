import math

# Step 1: Define the heights in inches.
height_with_milk_ft = 4
height_with_milk_in = 6
total_height_with_milk_inches = height_with_milk_ft * 12 + height_with_milk_in

height_without_milk_ft = 3
height_without_milk_in = 6
total_height_without_milk_inches = height_without_milk_ft * 12 + height_without_milk_in

# Step 2: Determine the proportion of the population for each height.
# The height of an offspring is determined by the father's genotype.
# If the father is 0/0, he cannot provide milk, and the offspring will be short.
# The frequency of the 0/0 genotype in the population is given as 0.5.
# Due to random mating, the probability of an individual's father being 0/0 is 0.5.
prop_short = 0.5

# The remaining proportion of fathers can provide milk, leading to tall offspring.
prop_tall = 1 - prop_short

# Step 3: Calculate the average height of the population.
# Average height = (proportion tall * tall height) + (proportion short * short height)
average_height = (prop_tall * total_height_with_milk_inches) + (prop_short * total_height_without_milk_inches)

# Step 4: Print the final equation and the result formatted to four significant figures.
# The result 48 needs to be formatted as 48.00 to have four significant figures.
print("The average adult height is the weighted average of the heights of individuals who did and did not receive milk.")
print(f"The equation is: ({prop_tall} * {total_height_with_milk_inches}) + ({prop_short} * {total_height_without_milk_inches})")
print(f"The average adult height of the population is {average_height:.2f} inches.")

# Final answer in the required format
# <<<48.00>>>