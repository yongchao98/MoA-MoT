# Number of each animal
num_horses = 15
num_humans = 11
num_lobsters = 5
num_shrimps = 11
num_ants = 6
num_sheeps = 8
num_scorpions = 1

# Number of legs for each animal
legs_per_horse = 4
legs_per_human = 2
legs_per_lobster = 10
legs_per_shrimp = 10
legs_per_ant = 6
legs_per_sheep = 4
legs_per_scorpion = 8

# Calculate total number of legs
total_legs = (num_horses * legs_per_horse +
              num_humans * legs_per_human +
              num_lobsters * legs_per_lobster +
              num_shrimps * legs_per_shrimp +
              num_ants * legs_per_ant +
              num_sheeps * legs_per_sheep +
              num_scorpions * legs_per_scorpion)

print(total_legs)