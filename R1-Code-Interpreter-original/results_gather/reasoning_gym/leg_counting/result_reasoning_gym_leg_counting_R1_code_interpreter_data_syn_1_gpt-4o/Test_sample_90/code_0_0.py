# Define the number of legs for each type of animal
legs_per_animal = {
    "starfish": 5,
    "leech": 0,
    "beetle": 6,
    "elephant": 4,
    "praying mantis": 6,
    "cat": 4
}

# Define the number of each type of animal
animal_counts = {
    "starfish": 2,
    "leech": 5,
    "beetle": 8,
    "elephant": 12,
    "praying mantis": 4,
    "cat": 9
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

# Output the result
print(total_legs)