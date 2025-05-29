# Number of animals
num_humans = 10
num_grasshoppers = 7
num_praying_mantises = 8
num_deers = 10

# Number of legs per animal
legs_per_human = 2
legs_per_grasshopper = 6
legs_per_praying_mantis = 6
legs_per_deer = 4

# Total number of legs
total_legs = (num_humans * legs_per_human) + \
             (num_grasshoppers * legs_per_grasshopper) + \
             (num_praying_mantises * legs_per_praying_mantis) + \
             (num_deers * legs_per_deer)

print(total_legs)