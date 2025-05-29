# Number of legs for each type of animal
legs_per_bird = 2
legs_per_insect = 6
legs_per_horse = 4
legs_per_elephant = 4
legs_per_praying_mantis = 6
legs_per_butterfly = 6

# Number of each type of animal
num_birds = 10
num_insects = 4
num_horses = 11
num_elephants = 7
num_praying_mantises = 9
num_butterflies = 14

# Calculate total number of legs
total_legs = (num_birds * legs_per_bird +
              num_insects * legs_per_insect +
              num_horses * legs_per_horse +
              num_elephants * legs_per_elephant +
              num_praying_mantises * legs_per_praying_mantis +
              num_butterflies * legs_per_butterfly)

print(total_legs)