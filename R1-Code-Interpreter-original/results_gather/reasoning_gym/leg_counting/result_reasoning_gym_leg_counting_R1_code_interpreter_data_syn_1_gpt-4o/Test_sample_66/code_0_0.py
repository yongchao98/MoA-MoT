# Define the number of legs for each type of animal
legs_per_animal = {
    "woodlouse": 14,
    "jellyfish": 0,
    "butterfly": 6,
    "grasshopper": 6,
    "flatworm": 0,
    "firefly": 6,
    "dog": 4,
    "beetle": 6
}

# Define the number of each type of animal
animal_counts = {
    "woodlouse": 14,
    "jellyfish": 11,
    "butterfly": 9,
    "grasshopper": 2,
    "flatworm": 13,
    "firefly": 8,
    "dog": 15,
    "beetle": 4
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

# Output the result
print(total_legs)