# Define the number of legs for each type of animal
legs_per_animal = {
    "human": 2,
    "scorpion": 8,
    "firefly": 6,
    "cow": 4,
    "praying_mantis": 6,
    "starfish": 5,
    "sheep": 4
}

# Define the number of each type of animal
animal_counts = {
    "human": 2,
    "scorpion": 8,
    "firefly": 7,
    "cow": 9,
    "praying_mantis": 9,
    "starfish": 4,
    "sheep": 8
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

# Output the result
print(total_legs)