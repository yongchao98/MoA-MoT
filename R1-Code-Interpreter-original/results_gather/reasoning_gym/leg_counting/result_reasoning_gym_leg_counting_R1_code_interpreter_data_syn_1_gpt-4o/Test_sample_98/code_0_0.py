# Number of animals
deers = 6
cockroaches = 11
praying_mantises = 9
sea_slugs = 15
butterflies = 11
giraffes = 12
sheeps = 2
cats = 9
humans = 11

# Number of legs per animal
legs_per_deer = 4
legs_per_cockroach = 6
legs_per_praying_mantis = 6
legs_per_sea_slug = 0
legs_per_butterfly = 6
legs_per_giraffe = 4
legs_per_sheep = 4
legs_per_cat = 4
legs_per_human = 2

# Total legs calculation
total_legs = (deers * legs_per_deer +
              cockroaches * legs_per_cockroach +
              praying_mantises * legs_per_praying_mantis +
              sea_slugs * legs_per_sea_slug +
              butterflies * legs_per_butterfly +
              giraffes * legs_per_giraffe +
              sheeps * legs_per_sheep +
              cats * legs_per_cat +
              humans * legs_per_human)

print(total_legs)