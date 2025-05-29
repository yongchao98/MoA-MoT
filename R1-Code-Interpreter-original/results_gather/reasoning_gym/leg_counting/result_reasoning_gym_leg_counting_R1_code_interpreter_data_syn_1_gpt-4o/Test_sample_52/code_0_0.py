# Define the number of legs for each type of animal
legs_per_animal = {
    "jellyfish": 0,
    "spider": 8,
    "snake": 0,
    "crab": 10,
    "insect": 6,
    "horse": 4,
    "elephant": 4
}

# Define the number of each type of animal
animal_counts = {
    "jellyfish": 8,
    "spider": 14,
    "snake": 9,
    "crab": 13,
    "insect": 12,
    "horse": 6,
    "elephant": 10
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

# Print the total number of legs
print(total_legs)