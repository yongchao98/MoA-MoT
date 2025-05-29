# Number of legs for each type of animal
legs_per_animal = {
    "praying_mantis": 6,
    "scorpion": 8,
    "firefly": 6,
    "crab": 10,
    "beetle": 6,
    "woodlouse": 14
}

# Number of each type of animal
animal_counts = {
    "praying_mantis": 12,
    "scorpion": 1,
    "firefly": 14,
    "crab": 15,
    "beetle": 1,
    "woodlouse": 15
}

# Calculate total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

print(total_legs)