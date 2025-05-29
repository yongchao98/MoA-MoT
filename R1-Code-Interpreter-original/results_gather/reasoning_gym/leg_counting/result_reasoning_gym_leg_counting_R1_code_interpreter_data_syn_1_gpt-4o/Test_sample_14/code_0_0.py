# Define the number of legs for each type of animal
legs_per_animal = {
    'wasp': 6,
    'jellyfish': 0,
    'elephant': 4,
    'crab': 10
}

# Define the number of each type of animal
animal_counts = {
    'wasp': 3,
    'jellyfish': 10,
    'elephant': 9,
    'crab': 13
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

# Output the result
print(total_legs)