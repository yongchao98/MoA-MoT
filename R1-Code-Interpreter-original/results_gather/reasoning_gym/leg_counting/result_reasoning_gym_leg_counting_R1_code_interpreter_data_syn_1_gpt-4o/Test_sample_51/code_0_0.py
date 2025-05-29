# Define the number of legs for each type of animal
legs_per_animal = {
    "praying_mantis": 6,
    "chicken": 2,
    "duck": 2,
    "bee": 6,
    "jellyfish": 0,
    "flatworm": 0,
    "firefly": 6,
    "insect": 6,
    "bird": 2,
    "cockroach": 6
}

# Define the number of each type of animal
animal_counts = {
    "praying_mantis": 3,
    "chicken": 5,
    "duck": 2,
    "bee": 14,
    "jellyfish": 1,
    "flatworm": 7,
    "firefly": 3,
    "insect": 6,
    "bird": 5,
    "cockroach": 9
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

# Output the result
print(total_legs)