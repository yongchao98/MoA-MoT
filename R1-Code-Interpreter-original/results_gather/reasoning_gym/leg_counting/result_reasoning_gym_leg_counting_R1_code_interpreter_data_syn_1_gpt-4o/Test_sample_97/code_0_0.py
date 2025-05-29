# Number of legs for each type of animal
legs_per_animal = {
    "jellyfish": 0,
    "wasp": 6,
    "elephant": 4,
    "ant": 6,
    "lobster": 10,
    "human": 2,
    "butterfly": 6
}

# Number of each type of animal
animal_counts = {
    "jellyfish": 9,
    "wasp": 3,
    "elephant": 3,
    "ant": 11,
    "lobster": 13,
    "human": 5,
    "butterfly": 10
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

print(total_legs)