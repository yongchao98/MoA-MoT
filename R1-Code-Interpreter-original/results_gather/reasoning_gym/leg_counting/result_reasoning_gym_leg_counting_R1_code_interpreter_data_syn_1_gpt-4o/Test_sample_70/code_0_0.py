# Number of legs for each type of animal
legs_per_butterfly = 6
legs_per_cat = 4
legs_per_shrimp = 10
legs_per_cow = 4
legs_per_bee = 6

# Number of each type of animal
num_butterflies = 5
num_cats = 10
num_shrimps = 7
num_cows = 10
num_bees = 10

# Calculate total number of legs
total_legs = (num_butterflies * legs_per_butterfly +
              num_cats * legs_per_cat +
              num_shrimps * legs_per_shrimp +
              num_cows * legs_per_cow +
              num_bees * legs_per_bee)

print(total_legs)