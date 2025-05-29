# Define the number of legs for each type of animal
legs_per_animal = {
    "cockroach": 6,
    "sheep": 4,
    "tiger": 4,
    "firefly": 6,
    "chicken": 2,
    "bee": 6,
    "beetle": 6,
    "jellyfish": 0,
    "praying mantis": 6,
    "ant": 6
}

# Define the number of each type of animal
animal_counts = {
    "cockroach": 4,
    "sheep": 5,
    "tiger": 2,
    "firefly": 14,
    "chicken": 11,
    "bee": 5,
    "beetle": 14,
    "jellyfish": 12,
    "praying mantis": 10,
    "ant": 1
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

# Output the result
print(total_legs)