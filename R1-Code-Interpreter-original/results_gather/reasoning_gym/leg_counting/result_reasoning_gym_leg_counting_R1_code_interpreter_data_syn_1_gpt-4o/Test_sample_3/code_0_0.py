# Number of each type of animal
num_ducks = 12
num_ants = 9
num_dogs = 13
num_sheep = 14
num_crickets = 15
num_tigers = 13
num_crabs = 5
num_sea_slugs = 11

# Number of legs for each type of animal
legs_duck = 2
legs_ant = 6
legs_dog = 4
legs_sheep = 4
legs_cricket = 6
legs_tiger = 4
legs_crab = 10
legs_sea_slug = 0

# Calculate total number of legs
total_legs = (num_ducks * legs_duck +
              num_ants * legs_ant +
              num_dogs * legs_dog +
              num_sheep * legs_sheep +
              num_crickets * legs_cricket +
              num_tigers * legs_tiger +
              num_crabs * legs_crab +
              num_sea_slugs * legs_sea_slug)

print(total_legs)