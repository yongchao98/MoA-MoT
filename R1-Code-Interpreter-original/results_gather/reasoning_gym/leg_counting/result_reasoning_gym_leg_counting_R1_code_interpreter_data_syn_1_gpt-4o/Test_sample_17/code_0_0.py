# Number of each type of animal
num_butterflies = 13
num_grasshoppers = 9
num_humans = 9
num_insects = 13

# Legs per animal
legs_per_butterfly = 6
legs_per_grasshopper = 6
legs_per_human = 2
legs_per_insect = 6

# Total legs calculation
total_legs = (num_butterflies * legs_per_butterfly) + \
             (num_grasshoppers * legs_per_grasshopper) + \
             (num_humans * legs_per_human) + \
             (num_insects * legs_per_insect)

print(total_legs)