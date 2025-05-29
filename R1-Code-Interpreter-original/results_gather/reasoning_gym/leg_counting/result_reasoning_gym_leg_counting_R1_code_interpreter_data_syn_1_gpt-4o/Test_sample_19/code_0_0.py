# Number of legs for each type of animal
legs_per_scorpion = 8
legs_per_shrimp = 10
legs_per_cat = 4
legs_per_giraffe = 4
legs_per_elephant = 4
legs_per_tiger = 4

# Number of each type of animal
num_scorpions = 15
num_shrimps = 14
num_cats = 7
num_giraffes = 9
num_elephants = 8
num_tigers = 12

# Calculate total number of legs
total_legs = (num_scorpions * legs_per_scorpion +
              num_shrimps * legs_per_shrimp +
              num_cats * legs_per_cat +
              num_giraffes * legs_per_giraffe +
              num_elephants * legs_per_elephant +
              num_tigers * legs_per_tiger)

print(total_legs)