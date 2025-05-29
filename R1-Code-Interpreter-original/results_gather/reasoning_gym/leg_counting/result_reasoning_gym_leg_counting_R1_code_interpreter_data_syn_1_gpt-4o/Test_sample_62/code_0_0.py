# Define the number of legs for each type of animal
legs_per_human = 2
legs_per_deer = 4
legs_per_starfish = 0
legs_per_cricket = 6
legs_per_sea_slug = 0

# Define the number of each type of animal
num_humans = 11
num_deers = 8
num_starfishs = 6
num_crickets = 5
num_sea_slugs = 7

# Calculate the total number of legs
total_legs = (num_humans * legs_per_human +
              num_deers * legs_per_deer +
              num_starfishs * legs_per_starfish +
              num_crickets * legs_per_cricket +
              num_sea_slugs * legs_per_sea_slug)

# Output the total number of legs
print(total_legs)