# Number of legs per animal
legs_per_cockroach = 6
legs_per_shrimp = 10
legs_per_starfish = 5
legs_per_firefly = 6
legs_per_spider = 8

# Number of each animal
num_cockroaches = 7
num_shrimps = 4
num_starfishes = 8
num_fireflies = 6
num_spiders = 10

# Calculate total number of legs
total_legs = (num_cockroaches * legs_per_cockroach +
              num_shrimps * legs_per_shrimp +
              num_starfishes * legs_per_starfish +
              num_fireflies * legs_per_firefly +
              num_spiders * legs_per_spider)

print(total_legs)