# Define the number of legs for each type of animal
legs_per_animal = {
    "cricket": 6,
    "praying_mantis": 6,
    "jellyfish": 0,
    "scorpion": 8,
    "firefly": 6,
    "butterfly": 6,
    "flatworm": 0,
    "cockroach": 6,
    "spider": 8,
    "cat": 4
}

# Define the number of each type of animal
animal_counts = {
    "cricket": 15,
    "praying_mantis": 2,
    "jellyfish": 13,
    "scorpion": 5,
    "firefly": 5,
    "butterfly": 8,
    "flatworm": 12,
    "cockroach": 14,
    "spider": 11,
    "cat": 6
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

# Output the total number of legs
print(total_legs)