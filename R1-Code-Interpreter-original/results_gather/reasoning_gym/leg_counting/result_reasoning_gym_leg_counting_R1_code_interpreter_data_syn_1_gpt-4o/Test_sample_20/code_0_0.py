# Number of legs for each type of animal
legs_per_animal = {
    "flatworm": 0,
    "scorpion": 8,
    "firefly": 6,
    "beetle": 6,
    "grasshopper": 6,
    "jellyfish": 0,
    "praying mantis": 6
}

# Number of each type of animal
animal_counts = {
    "flatworm": 3,
    "scorpion": 3,
    "firefly": 13,
    "beetle": 7,
    "grasshopper": 14,
    "jellyfish": 3,
    "praying mantis": 6
}

# Calculate total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

print(total_legs)