# Number of each type of animal
num_woodlouses = 11
num_praying_mantiss = 3
num_giraffes = 8

# Number of legs for each type of animal
legs_per_woodlouse = 14
legs_per_praying_mantis = 6
legs_per_giraffe = 4

# Calculate total number of legs
total_legs = (num_woodlouses * legs_per_woodlouse) + \
             (num_praying_mantiss * legs_per_praying_mantis) + \
             (num_giraffes * legs_per_giraffe)

print(total_legs)