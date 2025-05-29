# Define the number of legs for each type of animal
legs_per_animal = {
    "flatworms": 0,
    "sea slugs": 0,
    "sheep": 4,
    "chickens": 2,
    "snakes": 0,
    "cows": 4,
    "humans": 2,
    "ducks": 2
}

# Define the number of each type of animal
animal_counts = {
    "flatworms": 12,
    "sea slugs": 11,
    "sheep": 7,
    "chickens": 6,
    "snakes": 9,
    "cows": 5,
    "humans": 9,
    "ducks": 2
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

# Print the total number of legs
print(total_legs)