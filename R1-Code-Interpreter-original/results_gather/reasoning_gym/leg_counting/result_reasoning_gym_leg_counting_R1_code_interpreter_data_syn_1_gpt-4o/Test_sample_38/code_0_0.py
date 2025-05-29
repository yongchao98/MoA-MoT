# Number of animals
num_sheep = 6
num_praying_mantis = 2
num_cockroach = 3
num_duck = 14

# Legs per animal
legs_per_sheep = 4
legs_per_praying_mantis = 6
legs_per_cockroach = 6
legs_per_duck = 2

# Total legs calculation
total_legs = (num_sheep * legs_per_sheep) + \
             (num_praying_mantis * legs_per_praying_mantis) + \
             (num_cockroach * legs_per_cockroach) + \
             (num_duck * legs_per_duck)

print(total_legs)