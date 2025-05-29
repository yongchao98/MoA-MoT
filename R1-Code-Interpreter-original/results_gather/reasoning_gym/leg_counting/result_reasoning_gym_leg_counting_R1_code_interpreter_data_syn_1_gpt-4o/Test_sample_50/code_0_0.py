# Number of each type of animal
num_spiders = 11
num_scorpions = 5
num_sheep = 3
num_cows = 13
num_cats = 15

# Number of legs for each type of animal
legs_per_spider = 8
legs_per_scorpion = 8
legs_per_sheep = 4
legs_per_cow = 4
legs_per_cat = 4

# Calculate total number of legs
total_legs = (num_spiders * legs_per_spider) + \
             (num_scorpions * legs_per_scorpion) + \
             (num_sheep * legs_per_sheep) + \
             (num_cows * legs_per_cow) + \
             (num_cats * legs_per_cat)

print(total_legs)