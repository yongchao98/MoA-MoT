# Define the number of legs for each type of animal
legs_per_animal = {
    "cow": 4,
    "shrimp": 10,
    "lion": 4,
    "insect": 6,
    "sheep": 4,
    "crab": 10,
    "wasp": 6,
    "butterfly": 6,
    "snake": 0,
    "giraffe": 4
}

# Define the number of each type of animal
animal_counts = {
    "cow": 6,
    "shrimp": 6,
    "lion": 2,
    "insect": 7,
    "sheep": 8,
    "crab": 6,
    "wasp": 3,
    "butterfly": 6,
    "snake": 11,
    "giraffe": 1
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

# Output the total number of legs
print(total_legs)