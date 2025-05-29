# Define the number of legs for each type of animal
legs_per_animal = {
    "birds": 2,
    "snakes": 0,
    "starfish": 0,
    "crabs": 10,
    "ants": 6,
    "sheep": 4,
    "chickens": 2,
    "crickets": 6
}

# Define the number of each type of animal
animal_counts = {
    "birds": 6,
    "snakes": 3,
    "starfish": 13,
    "crabs": 15,
    "ants": 3,
    "sheep": 2,
    "chickens": 11,
    "crickets": 7
}

# Calculate the total number of legs
total_legs = sum(animal_counts[animal] * legs_per_animal[animal] for animal in animal_counts)

# Output the result
print(total_legs)