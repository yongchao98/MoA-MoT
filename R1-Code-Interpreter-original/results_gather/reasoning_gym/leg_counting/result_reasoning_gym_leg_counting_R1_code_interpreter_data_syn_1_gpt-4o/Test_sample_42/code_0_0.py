# Define the number of legs for each type of animal
legs_per_insect = 6
legs_per_human = 2
legs_per_starfish = 0
legs_per_cow = 4

# Define the number of each type of animal
num_insects = 15
num_humans = 14
num_starfish = 15
num_cows = 4

# Calculate the total number of legs
total_legs = (num_insects * legs_per_insect) + (num_humans * legs_per_human) + (num_starfish * legs_per_starfish) + (num_cows * legs_per_cow)

# Output the total number of legs
print(total_legs)