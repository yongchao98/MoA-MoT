# Define the number of legs for each type of animal
legs_per_animal = {
    "sea_slug": 0,
    "elephant": 4,
    "woodlouse": 14,
    "bee": 6,
    "beetle": 6,
    "leech": 0,
    "lion": 4,
    "dog": 4
}

# Define the number of each animal
animal_counts = {
    "sea_slug": 6,
    "elephant": 2,
    "woodlouse": 10,
    "bee": 9,
    "beetle": 14,
    "leech": 3,
    "lion": 2,
    "dog": 5
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

# Output the total number of legs
print(total_legs)