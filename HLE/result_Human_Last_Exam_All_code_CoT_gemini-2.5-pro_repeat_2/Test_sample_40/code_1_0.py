import math

# Step 1: Define the given parameters.
# The frequency of the 0/0 genotype in the population is 0.5.
# This is the frequency of fathers who cannot provide milk.
freq_no_milk_provider = 0.5

# The frequency of fathers who can provide milk (genotypes 1/1 or 1/0).
freq_milk_provider = 1 - freq_no_milk_provider

# Height for individuals with milk during upbringing (4 feet 6 inches).
height_with_milk_inches = 4 * 12 + 6

# Height for individuals without milk during upbringing (3 feet 6 inches).
height_without_milk_inches = 3 * 12 + 6

# Step 2: Calculate the average height of the population.
# This is a weighted average based on the proportions of the two height groups.
# The proportion of the population in each height group is determined by the frequency
# of their father's genotype.
proportion_tall = freq_milk_provider
proportion_short = freq_no_milk_provider

average_height = (proportion_tall * height_with_milk_inches) + \
                 (proportion_short * height_without_milk_inches)

# Step 3: Print the explanation and the final equation.
print("To find the average height of the population, we calculate a weighted average.")
print("The population is divided into two groups based on the milk provided during their upbringing, which depends on their father's genotype.")
print(f"\n1. Individuals with milk (tall group):")
print(f"   - Height: 4'6\" = {height_with_milk_inches} inches")
print(f"   - Proportion of population: This corresponds to the frequency of fathers who can provide milk, which is {proportion_tall:.1f}.")

print(f"\n2. Individuals without milk (short group):")
print(f"   - Height: 3'6\" = {height_without_milk_inches} inches")
print(f"   - Proportion of population: This corresponds to the frequency of fathers with the 0/0 genotype, which is {proportion_short:.1f}.")

print("\n3. The average height is calculated using the following equation:")
# The user wants to see the final equation with the numbers.
print(f"Average Height = ({proportion_tall} * {height_with_milk_inches}) + ({proportion_short} * {height_without_milk_inches})")

# Step 4: Display the final result formatted to four significant figures.
# 48.00 has four significant figures.
formatted_average_height = f"{average_height:.2f}"
print(f"\nThe calculated average adult height is {formatted_average_height} inches.")

# The final answer in the requested format.
print(f"\n<<<{formatted_average_height}>>>")