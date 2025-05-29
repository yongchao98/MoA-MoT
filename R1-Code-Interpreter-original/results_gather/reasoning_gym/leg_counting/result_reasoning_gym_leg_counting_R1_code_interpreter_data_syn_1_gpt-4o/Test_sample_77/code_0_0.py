# Number of animals
elephants = 8
praying_mantises = 5
cockroaches = 10
flatworms = 15
birds = 9
woodlice = 2
tigers = 13

# Number of legs per animal
legs_per_elephant = 4
legs_per_praying_mantis = 6
legs_per_cockroach = 6
legs_per_flatworm = 0
legs_per_bird = 2
legs_per_woodlouse = 14
legs_per_tiger = 4

# Total number of legs
total_legs = (elephants * legs_per_elephant +
              praying_mantises * legs_per_praying_mantis +
              cockroaches * legs_per_cockroach +
              flatworms * legs_per_flatworm +
              birds * legs_per_bird +
              woodlice * legs_per_woodlouse +
              tigers * legs_per_tiger)

print(total_legs)