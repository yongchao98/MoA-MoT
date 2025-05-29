# Number of animals
bees = 12
horses = 7
cows = 9
elephants = 11
giraffes = 12
ducks = 9
woodlouses = 2
jellyfishs = 10
spiders = 8

# Number of legs per animal
legs_per_bee = 6
legs_per_horse = 4
legs_per_cow = 4
legs_per_elephant = 4
legs_per_giraffe = 4
legs_per_duck = 2
legs_per_woodlouse = 14
legs_per_jellyfish = 0
legs_per_spider = 8

# Total number of legs
total_legs = (bees * legs_per_bee +
              horses * legs_per_horse +
              cows * legs_per_cow +
              elephants * legs_per_elephant +
              giraffes * legs_per_giraffe +
              ducks * legs_per_duck +
              woodlouses * legs_per_woodlouse +
              jellyfishs * legs_per_jellyfish +
              spiders * legs_per_spider)

print(total_legs)