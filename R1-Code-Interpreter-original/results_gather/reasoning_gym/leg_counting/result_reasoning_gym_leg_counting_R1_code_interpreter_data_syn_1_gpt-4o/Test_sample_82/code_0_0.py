# Define the number of legs for each type of animal
legs_per_animal = {
    "woodlouse": 14,
    "bird": 2,
    "ant": 6,
    "giraffe": 4,
    "flatworm": 0,
    "dog": 4
}

# Define the number of each type of animal
animal_counts = {
    "woodlouse": 9,
    "bird": 9,
    "ant": 11,
    "giraffe": 9,
    "flatworm": 5,
    "dog": 10
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

# Output the result
print(total_legs)