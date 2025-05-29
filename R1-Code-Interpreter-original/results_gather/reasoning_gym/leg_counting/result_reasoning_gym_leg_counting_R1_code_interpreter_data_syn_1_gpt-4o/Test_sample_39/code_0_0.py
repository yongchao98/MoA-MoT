# Define the number of legs for each type of animal
legs_per_animal = {
    "scorpion": 8,
    "snake": 0,
    "tiger": 4,
    "dog": 4,
    "horse": 4,
    "lobster": 10,
    "jellyfish": 0,
    "bird": 2
}

# Define the number of each type of animal
animal_counts = {
    "scorpion": 9,
    "snake": 1,
    "tiger": 10,
    "dog": 4,
    "horse": 9,
    "lobster": 1,
    "jellyfish": 14,
    "bird": 12
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

# Output the result
print(total_legs)