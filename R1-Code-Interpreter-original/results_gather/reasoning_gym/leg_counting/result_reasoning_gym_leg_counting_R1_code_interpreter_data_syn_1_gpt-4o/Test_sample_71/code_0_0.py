# Define the number of legs for each type of animal
legs_per_animal = {
    "lion": 4,
    "leech": 0,
    "human": 2,
    "bee": 6,
    "praying_mantis": 6,
    "ant": 6,
    "cricket": 6,
    "spider": 8,
    "flatworm": 0,
    "giraffe": 4
}

# Define the number of each type of animal
animal_counts = {
    "lion": 2,
    "leech": 3,
    "human": 4,
    "bee": 10,
    "praying_mantis": 14,
    "ant": 1,
    "cricket": 1,
    "spider": 4,
    "flatworm": 2,
    "giraffe": 12
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

# Output the total number of legs
print(total_legs)