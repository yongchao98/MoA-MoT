# Define the number of legs for each type of animal
legs_per_animal = {
    "cricket": 6,
    "bird": 2,
    "leech": 0,
    "beetle": 6,
    "praying mantis": 6,
    "crab": 10,
    "sea slug": 0,
    "starfish": 0,
    "ant": 6
}

# Define the number of each type of animal
animal_counts = {
    "cricket": 3,
    "bird": 10,
    "leech": 6,
    "beetle": 6,
    "praying mantis": 5,
    "crab": 12,
    "sea slug": 6,
    "starfish": 10,
    "ant": 12
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

# Output the total number of legs
print(total_legs)