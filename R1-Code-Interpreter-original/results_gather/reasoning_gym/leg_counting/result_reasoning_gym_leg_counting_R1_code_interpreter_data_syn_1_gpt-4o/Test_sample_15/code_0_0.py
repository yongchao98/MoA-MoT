# Define the number of legs for each type of animal
legs_per_animal = {
    "leech": 0,
    "cricket": 6,
    "praying_mantis": 6,
    "ant": 6,
    "tiger": 4,
    "bee": 6,
    "horse": 4,
    "starfish": 0
}

# Define the number of each type of animal
animal_counts = {
    "leech": 1,
    "cricket": 13,
    "praying_mantis": 7,
    "ant": 11,
    "tiger": 1,
    "bee": 2,
    "horse": 1,
    "starfish": 5
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

# Output the result
print(total_legs)