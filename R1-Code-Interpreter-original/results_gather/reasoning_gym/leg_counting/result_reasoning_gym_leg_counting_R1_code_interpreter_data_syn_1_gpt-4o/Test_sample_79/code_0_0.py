# Define the number of legs for each type of animal
legs_per_insect = 6
legs_per_spider = 8
legs_per_horse = 4
legs_per_wasp = 6
legs_per_grasshopper = 6
legs_per_cat = 4

# Define the number of each type of animal
num_insects = 11
num_spiders = 4
num_horses = 12
num_wasps = 7
num_grasshoppers = 9
num_cats = 4

# Calculate the total number of legs
total_legs = (num_insects * legs_per_insect +
              num_spiders * legs_per_spider +
              num_horses * legs_per_horse +
              num_wasps * legs_per_wasp +
              num_grasshoppers * legs_per_grasshopper +
              num_cats * legs_per_cat)

# Output the total number of legs
print(total_legs)