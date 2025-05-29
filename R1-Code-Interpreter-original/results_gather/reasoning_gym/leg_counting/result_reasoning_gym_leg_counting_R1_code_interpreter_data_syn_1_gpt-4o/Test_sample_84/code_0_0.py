# Define the number of legs for each type of animal
legs_per_animal = {
    'snake': 0,
    'human': 2,
    'duck': 2,
    'firefly': 6,
    'crab': 10,
    'giraffe': 4,
    'lobster': 10
}

# Define the number of each type of animal
animal_counts = {
    'snake': 7,
    'human': 15,
    'duck': 10,
    'firefly': 7,
    'crab': 14,
    'giraffe': 15,
    'lobster': 10
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

# Output the result
print(total_legs)