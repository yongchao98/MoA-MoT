# Define the number of legs for each type of animal
legs_per_animal = {
    "cricket": 6,
    "firefly": 6,
    "starfish": 5,
    "duck": 2,
    "woodlouse": 14,
    "human": 2
}

# Define the number of each type of animal
animal_counts = {
    "cricket": 5,
    "firefly": 5,
    "starfish": 10,
    "duck": 1,
    "woodlouse": 12,
    "human": 11
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

# Output the result
print(total_legs)