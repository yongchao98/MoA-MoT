# Define the number of legs for each type of animal
legs_per_scorpion = 8
legs_per_jellyfish = 0
legs_per_beetle = 6
legs_per_cockroach = 6

# Define the number of each type of animal
num_scorpions = 9
num_jellyfishs = 13
num_beetles = 2
num_cockroachs = 11

# Calculate the total number of legs
total_legs = (num_scorpions * legs_per_scorpion +
              num_jellyfishs * legs_per_jellyfish +
              num_beetles * legs_per_beetle +
              num_cockroachs * legs_per_cockroach)

# Output the total number of legs
print(total_legs)