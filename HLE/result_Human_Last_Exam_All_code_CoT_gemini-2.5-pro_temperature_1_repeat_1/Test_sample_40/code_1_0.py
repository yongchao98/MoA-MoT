import math

# Step 1: Define the given information
# The problem states that half the population has the 0/0 genotype.
# Milk availability depends on the father's genotype. Due to random mating,
# the frequency of fathers with a 0/0 genotype is the same as the population frequency.
prop_without_milk = 0.5
prop_with_milk = 1 - prop_without_milk

# Heights in feet and inches
height_with_milk_ft, height_with_milk_in = 4, 6
height_without_milk_ft, height_without_milk_in = 3, 6

# Step 2: Convert heights to inches
height_with_milk_inches = height_with_milk_ft * 12 + height_with_milk_in
height_without_milk_inches = height_without_milk_ft * 12 + height_without_milk_in

# Step 3: Calculate the average height of the population
# This is a weighted average of the two possible heights.
average_height = (prop_with_milk * height_with_milk_inches) + (prop_without_milk * height_without_milk_inches)

# Step 4: Print the equation and the final answer
# The final answer is rounded to four significant figures.
print("The average adult height is calculated as a weighted average based on the proportion of the population that receives milk during development.")
print("The proportion of the population that gets milk is equal to the proportion of fathers who can digest milk.")
print(f"Proportion with milk: {prop_with_milk:.4f}")
print(f"Proportion without milk: {prop_without_milk:.4f}")
print(f"Height with milk: {height_with_milk_inches} inches")
print(f"Height without milk: {height_without_milk_inches} inches")
print("\nFinal Calculation:")
# Displaying the equation with each number as requested
print(f"({prop_with_milk:.4f} * {height_with_milk_inches}) + ({prop_without_milk:.4f} * {height_without_milk_inches}) = {average_height:.4f}")
print(f"\nThe average adult height of the population is {average_height:.4f} inches.")

# Final answer in the required format
print(f"\n<<<{average_height:.4f}>>>")