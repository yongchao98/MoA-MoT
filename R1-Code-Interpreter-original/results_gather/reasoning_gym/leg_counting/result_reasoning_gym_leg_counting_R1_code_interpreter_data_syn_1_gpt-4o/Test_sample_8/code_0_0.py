# Define the number of legs for each animal
legs_per_animal = {
    "crab": 10,
    "jellyfish": 0,
    "giraffe": 4,
    "sheep": 4,
    "lion": 4,
    "duck": 2,
    "bee": 6,
    "cow": 4,
    "spider": 8,
    "beetle": 6
}

# Define the number of each animal
animal_counts = {
    "crab": 8,
    "jellyfish": 14,
    "giraffe": 1,
    "sheep": 9,
    "lion": 15,
    "duck": 12,
    "bee": 3,
    "cow": 1,
    "spider": 11,
    "beetle": 4
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

# Output the total number of legs
print(total_legs)