# Number of legs for each type of animal
legs_per_mantis = 6
legs_per_cow = 4
legs_per_wasp = 6
legs_per_butterfly = 6
legs_per_cricket = 6

# Number of each type of animal
num_mantises = 9
num_cows = 13
num_wasps = 2
num_butterflies = 9
num_crickets = 2

# Calculate total number of legs
total_legs = (num_mantises * legs_per_mantis +
              num_cows * legs_per_cow +
              num_wasps * legs_per_wasp +
              num_butterflies * legs_per_butterfly +
              num_crickets * legs_per_cricket)

print(total_legs)