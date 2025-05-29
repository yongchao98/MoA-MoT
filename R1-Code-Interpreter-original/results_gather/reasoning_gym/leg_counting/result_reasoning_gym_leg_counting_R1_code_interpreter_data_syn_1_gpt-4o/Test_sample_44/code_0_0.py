# Define the number of legs for each type of animal
legs_per_animal = {
    "praying_mantis": 6,
    "bee": 6,
    "horse": 4,
    "starfish": 5,
    "flatworm": 0,
    "grasshopper": 6,
    "spider": 8,
    "human": 2
}

# Define the number of each type of animal
animal_counts = {
    "praying_mantis": 5,
    "bee": 5,
    "horse": 13,
    "starfish": 14,
    "flatworm": 8,
    "grasshopper": 8,
    "spider": 4,
    "human": 2
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

# Output the result
print(total_legs)