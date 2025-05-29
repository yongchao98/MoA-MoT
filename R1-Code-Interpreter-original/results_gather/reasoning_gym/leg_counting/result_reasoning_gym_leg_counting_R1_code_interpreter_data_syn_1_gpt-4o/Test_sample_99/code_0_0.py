# Number of animals
beetles = 15
scorpions = 2
chickens = 6
wasps = 10
sheep = 10
sea_slugs = 2
horses = 12
bees = 3

# Number of legs per animal
legs_per_beetle = 6
legs_per_scorpion = 8
legs_per_chicken = 2
legs_per_wasp = 6
legs_per_sheep = 4
legs_per_sea_slug = 0
legs_per_horse = 4
legs_per_bee = 6

# Total number of legs
total_legs = (beetles * legs_per_beetle +
              scorpions * legs_per_scorpion +
              chickens * legs_per_chicken +
              wasps * legs_per_wasp +
              sheep * legs_per_sheep +
              sea_slugs * legs_per_sea_slug +
              horses * legs_per_horse +
              bees * legs_per_bee)

print(total_legs)