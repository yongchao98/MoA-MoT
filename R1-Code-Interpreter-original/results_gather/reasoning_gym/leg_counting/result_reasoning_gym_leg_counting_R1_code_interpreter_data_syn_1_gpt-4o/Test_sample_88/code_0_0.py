# Number of animals
horses = 14
praying_mantises = 5
chickens = 12
insects = 6
lions = 6
elephants = 12

# Number of legs per animal
legs_per_horse = 4
legs_per_praying_mantis = 6
legs_per_chicken = 2
legs_per_insect = 6
legs_per_lion = 4
legs_per_elephant = 4

# Total legs calculation
total_legs = (horses * legs_per_horse +
              praying_mantises * legs_per_praying_mantis +
              chickens * legs_per_chicken +
              insects * legs_per_insect +
              lions * legs_per_lion +
              elephants * legs_per_elephant)

print(total_legs)