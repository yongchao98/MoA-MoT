# Define the number of legs for each type of animal
legs_per_animal = {
    "woodlouse": 14,
    "butterfly": 6,
    "human": 2,
    "chicken": 2,
    "shrimp": 10,
    "dog": 4,
    "crab": 10,
    "horse": 4,
    "flatworm": 0,
    "wasp": 6
}

# Define the number of each type of animal
animal_counts = {
    "woodlouse": 3,
    "butterfly": 10,
    "human": 1,
    "chicken": 10,
    "shrimp": 9,
    "dog": 15,
    "crab": 2,
    "horse": 11,
    "flatworm": 10,
    "wasp": 14
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

# Output the total number of legs
print(total_legs)