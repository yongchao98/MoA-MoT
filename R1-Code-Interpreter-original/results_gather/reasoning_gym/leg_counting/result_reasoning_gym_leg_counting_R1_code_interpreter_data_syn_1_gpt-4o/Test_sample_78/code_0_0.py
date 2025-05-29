# Define the number of legs for each type of animal
legs_per_animal = {
    "giraffe": 4,
    "grasshopper": 6,
    "firefly": 6,
    "deer": 4,
    "leech": 0,
    "shrimp": 10,
    "tiger": 4,
    "cat": 4,
    "beetle": 6
}

# Define the number of each type of animal
animal_counts = {
    "giraffe": 1,
    "grasshopper": 2,
    "firefly": 14,
    "deer": 9,
    "leech": 15,
    "shrimp": 2,
    "tiger": 3,
    "cat": 8,
    "beetle": 11
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

# Output the total number of legs
print(total_legs)