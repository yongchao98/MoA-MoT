# Plan: Calculate the average population height based on the father's genotype.
# The father's ability to provide milk determines the offspring's height.

# Step 1: Define the proportions of fathers with and without milk-digesting ability.
# We are told that half the population has the 0/0 genotype.
# Since mating is random, the probability that a child has a father with the 0/0
# genotype is the same as the frequency of that genotype in the population.
# These fathers cannot provide milk.
proportion_no_milk = 0.5

# The rest of the fathers can provide milk.
proportion_with_milk = 1 - proportion_no_milk

# Step 2: Define the heights in inches.
# Height without milk is 3 feet 6 inches.
height_no_milk_inches = 3 * 12 + 6

# Height with milk is 1 foot (12 inches) taller.
height_with_milk_inches = height_no_milk_inches + 12

# Step 3: Calculate the weighted average height of the population.
# Average Height = (Proportion of population raised without milk * Height without milk) +
#                  (Proportion of population raised with milk * Height with milk)
average_height = (proportion_no_milk * height_no_milk_inches) + (proportion_with_milk * height_with_milk_inches)

# Step 4: Print the final equation and the result formatted to four significant figures.
# Using integer values for the printout to match the input numbers.
print("The final calculation is: ({:.1f} * {:d}) + ({:.1f} * {:d})".format(proportion_no_milk, height_no_milk_inches, proportion_with_milk, height_with_milk_inches))
print("The average adult height of the population is {:.2f} inches.".format(average_height))
