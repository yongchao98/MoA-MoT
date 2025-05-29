# Number of legs for each type of animal
legs_per_cricket = 6
legs_per_beetle = 6
legs_per_flatworm = 0
legs_per_sea_slug = 0
legs_per_jellyfish = 0
legs_per_crab = 10

# Quantities of each type of animal
num_crickets = 15
num_beetles = 12
num_flatworms = 14
num_sea_slugs = 5
num_jellyfish = 6
num_crabs = 1

# Calculate total number of legs
total_legs = (num_crickets * legs_per_cricket +
              num_beetles * legs_per_beetle +
              num_flatworms * legs_per_flatworm +
              num_sea_slugs * legs_per_sea_slug +
              num_jellyfish * legs_per_jellyfish +
              num_crabs * legs_per_crab)

print(total_legs)