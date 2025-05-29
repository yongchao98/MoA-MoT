# Define the number of legs for each type of animal
legs_per_tiger = 4
legs_per_duck = 2
legs_per_butterfly = 6
legs_per_leech = 0

# Define the number of each type of animal
num_tigers = 3
num_ducks = 14
num_butterflies = 7
num_leeches = 6

# Calculate the total number of legs
total_legs = (num_tigers * legs_per_tiger) + (num_ducks * legs_per_duck) + \
             (num_butterflies * legs_per_butterfly) + (num_leeches * legs_per_leech)

# Print the total number of legs
print(total_legs)