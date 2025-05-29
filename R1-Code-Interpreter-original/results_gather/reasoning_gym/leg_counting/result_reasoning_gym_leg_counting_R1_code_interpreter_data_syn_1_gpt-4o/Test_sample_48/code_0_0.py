# Number of legs for each type of animal
legs_per_animal = {
    "wasp": 6,
    "lobster": 10,
    "grasshopper": 6,
    "beetle": 6,
    "crab": 10,
    "woodlouse": 14,
    "cow": 4
}

# Number of each type of animal
animal_counts = {
    "wasp": 1,
    "lobster": 9,
    "grasshopper": 13,
    "beetle": 14,
    "crab": 7,
    "woodlouse": 15,
    "cow": 4
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

print(total_legs)