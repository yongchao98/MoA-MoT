# Number of legs for each type of animal
legs_per_animal = {
    "shrimp": 10,
    "snake": 0,
    "cow": 4,
    "flatworm": 0,
    "grasshopper": 6,
    "deer": 4,
    "wasp": 6
}

# Number of each type of animal
animal_counts = {
    "shrimp": 6,
    "snake": 1,
    "cow": 6,
    "flatworm": 3,
    "grasshopper": 7,
    "deer": 5,
    "wasp": 14
}

# Calculate total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

print(total_legs)