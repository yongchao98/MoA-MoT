# Number of each type of animal
num_leeches = 3
num_giraffes = 9
num_praying_mantises = 3
num_starfish = 3
num_butterflies = 1

# Number of legs per type of animal
legs_per_leech = 0
legs_per_giraffe = 4
legs_per_praying_mantis = 6
legs_per_starfish = 0
legs_per_butterfly = 6

# Total number of legs
total_legs = (num_leeches * legs_per_leech +
              num_giraffes * legs_per_giraffe +
              num_praying_mantises * legs_per_praying_mantis +
              num_starfish * legs_per_starfish +
              num_butterflies * legs_per_butterfly)

print(total_legs)