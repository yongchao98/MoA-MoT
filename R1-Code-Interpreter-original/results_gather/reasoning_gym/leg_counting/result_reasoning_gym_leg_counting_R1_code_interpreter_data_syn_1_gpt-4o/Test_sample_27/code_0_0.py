# Number of animals
chickens = 11
praying_mantises = 10
sea_slugs = 1
ducks = 15
starfish = 14
lobsters = 8

# Number of legs per animal
legs_per_chicken = 2
legs_per_praying_mantis = 6
legs_per_sea_slug = 0
legs_per_duck = 2
legs_per_starfish = 0
legs_per_lobster = 10

# Total number of legs
total_legs = (chickens * legs_per_chicken +
              praying_mantises * legs_per_praying_mantis +
              sea_slugs * legs_per_sea_slug +
              ducks * legs_per_duck +
              starfish * legs_per_starfish +
              lobsters * legs_per_lobster)

print(total_legs)