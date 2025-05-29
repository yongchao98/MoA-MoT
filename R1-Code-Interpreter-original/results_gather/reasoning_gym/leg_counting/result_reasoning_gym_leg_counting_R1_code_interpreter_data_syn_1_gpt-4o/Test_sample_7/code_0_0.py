# Number of each animal
num_tigers = 12
num_bees = 10
num_sea_slugs = 1
num_horses = 11
num_elephants = 13
num_wasps = 5

# Number of legs for each type of animal
legs_per_tiger = 4
legs_per_bee = 6
legs_per_sea_slug = 0
legs_per_horse = 4
legs_per_elephant = 4
legs_per_wasp = 6

# Calculate total number of legs
total_legs = (num_tigers * legs_per_tiger +
              num_bees * legs_per_bee +
              num_sea_slugs * legs_per_sea_slug +
              num_horses * legs_per_horse +
              num_elephants * legs_per_elephant +
              num_wasps * legs_per_wasp)

print(total_legs)