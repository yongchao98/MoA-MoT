# Number of animals
cows = 4
lobsters = 3
crickets = 5
birds = 13
chickens = 5
shrimps = 4
insects = 8

# Number of legs per animal
legs_per_cow = 4
legs_per_lobster = 10
legs_per_cricket = 6
legs_per_bird = 2
legs_per_chicken = 2
legs_per_shrimp = 10
legs_per_insect = 6

# Total number of legs
total_legs = (cows * legs_per_cow +
              lobsters * legs_per_lobster +
              crickets * legs_per_cricket +
              birds * legs_per_bird +
              chickens * legs_per_chicken +
              shrimps * legs_per_shrimp +
              insects * legs_per_insect)

print(total_legs)