# Number of animals
num_sheep = 15
num_woodlouse = 6
num_scorpion = 14
num_elephant = 11
num_cockroach = 7
num_praying_mantis = 15
num_tiger = 4
num_crab = 2
num_deer = 15

# Number of legs per animal
legs_sheep = 4
legs_woodlouse = 14
legs_scorpion = 8
legs_elephant = 4
legs_cockroach = 6
legs_praying_mantis = 6
legs_tiger = 4
legs_crab = 10
legs_deer = 4

# Total legs calculation
total_legs = (num_sheep * legs_sheep +
              num_woodlouse * legs_woodlouse +
              num_scorpion * legs_scorpion +
              num_elephant * legs_elephant +
              num_cockroach * legs_cockroach +
              num_praying_mantis * legs_praying_mantis +
              num_tiger * legs_tiger +
              num_crab * legs_crab +
              num_deer * legs_deer)

print(total_legs)