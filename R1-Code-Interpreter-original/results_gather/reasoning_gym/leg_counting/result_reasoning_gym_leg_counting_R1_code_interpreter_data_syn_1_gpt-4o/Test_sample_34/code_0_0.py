# Define the number of legs for each type of animal
legs_per_starfish = 5
legs_per_cricket = 6
legs_per_giraffe = 4
legs_per_cockroach = 6

# Define the number of each type of animal
num_starfish = 5
num_crickets = 5
num_giraffes = 14
num_cockroaches = 13

# Calculate the total number of legs
total_legs = (num_starfish * legs_per_starfish +
              num_crickets * legs_per_cricket +
              num_giraffes * legs_per_giraffe +
              num_cockroaches * legs_per_cockroach)

# Output the total number of legs
print(total_legs)