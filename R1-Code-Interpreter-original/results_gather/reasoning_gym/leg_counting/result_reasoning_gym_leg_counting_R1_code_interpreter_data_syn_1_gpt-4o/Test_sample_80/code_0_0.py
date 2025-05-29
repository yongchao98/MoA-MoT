# Number of legs for each type of animal
legs_per_animal = {
    "starfish": 5,
    "praying_mantis": 6,
    "grasshopper": 6,
    "cat": 4,
    "beetle": 6,
    "elephant": 4,
    "tiger": 4,
    "spider": 8,
    "snake": 0,
    "horse": 4
}

# Number of each type of animal
animal_counts = {
    "starfish": 10,
    "praying_mantis": 4,
    "grasshopper": 12,
    "cat": 8,
    "beetle": 8,
    "elephant": 5,
    "tiger": 6,
    "spider": 4,
    "snake": 4,
    "horse": 12
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

print(total_legs)