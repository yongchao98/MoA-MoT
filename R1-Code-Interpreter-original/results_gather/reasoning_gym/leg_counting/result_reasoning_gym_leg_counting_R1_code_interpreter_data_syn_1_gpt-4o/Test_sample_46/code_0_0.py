# Number of each animal
birds = 15
starfish = 11
tigers = 7
crickets = 13
jellyfish = 11
dogs = 4

# Legs per animal
legs_per_bird = 2
legs_per_starfish = 0
legs_per_tiger = 4
legs_per_cricket = 6
legs_per_jellyfish = 0
legs_per_dog = 4

# Total legs calculation
total_legs = (birds * legs_per_bird) + (starfish * legs_per_starfish) + \
             (tigers * legs_per_tiger) + (crickets * legs_per_cricket) + \
             (jellyfish * legs_per_jellyfish) + (dogs * legs_per_dog)

print(total_legs)