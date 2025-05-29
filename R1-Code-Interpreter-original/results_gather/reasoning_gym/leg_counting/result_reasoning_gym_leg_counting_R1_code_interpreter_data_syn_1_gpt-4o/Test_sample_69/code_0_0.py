# Number of legs for each type of animal
legs_per_animal = {
    "lobster": 10,
    "butterfly": 6,
    "praying_mantis": 6,
    "flatworm": 0,
    "scorpion": 8,
    "lion": 4
}

# Number of each type of animal
animal_counts = {
    "lobster": 10,
    "butterfly": 4,
    "praying_mantis": 15,
    "flatworm": 5,
    "scorpion": 8,
    "lion": 12
}

# Calculate total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

print(total_legs)