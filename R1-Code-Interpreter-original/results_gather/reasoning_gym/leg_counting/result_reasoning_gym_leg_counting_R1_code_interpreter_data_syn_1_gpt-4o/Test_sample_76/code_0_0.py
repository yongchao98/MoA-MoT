# Define the number of legs for each type of animal
legs_per_animal = {
    'chicken': 2,
    'bee': 6,
    'deer': 4,
    'horse': 4,
    'starfish': 0
}

# Define the number of each type of animal
animal_counts = {
    'chicken': 4,
    'bee': 3,
    'deer': 4,
    'horse': 12,
    'starfish': 14
}

# Calculate the total number of legs
total_legs = sum(animal_counts[animal] * legs_per_animal[animal] for animal in animal_counts)

# Output the result
print(total_legs)