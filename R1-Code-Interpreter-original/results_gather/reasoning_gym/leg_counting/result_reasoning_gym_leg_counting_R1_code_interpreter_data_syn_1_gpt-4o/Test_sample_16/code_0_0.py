# Number of legs for each type of animal
legs_per_animal = {
    "duck": 2,
    "shrimp": 10,
    "lobster": 10,
    "starfish": 0,
    "scorpion": 8,
    "horse": 4,
    "butterfly": 6,
    "bee": 6
}

# Number of each type of animal
animal_counts = {
    "duck": 7,
    "shrimp": 3,
    "lobster": 3,
    "starfish": 13,
    "scorpion": 14,
    "horse": 11,
    "butterfly": 2,
    "bee": 8
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

print(total_legs)