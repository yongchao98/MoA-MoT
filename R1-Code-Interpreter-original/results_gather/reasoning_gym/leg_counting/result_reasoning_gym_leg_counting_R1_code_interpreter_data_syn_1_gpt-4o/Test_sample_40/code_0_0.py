# Define the number of legs for each type of animal
legs_per_cow = 4
legs_per_butterfly = 6
legs_per_lobster = 10
legs_per_beetle = 6

# Define the number of each type of animal
num_cows = 5
num_butterflies = 12
num_lobsters = 12
num_beetles = 11

# Calculate the total number of legs
total_legs = (num_cows * legs_per_cow) + \
             (num_butterflies * legs_per_butterfly) + \
             (num_lobsters * legs_per_lobster) + \
             (num_beetles * legs_per_beetle)

# Output the total number of legs
print(total_legs)