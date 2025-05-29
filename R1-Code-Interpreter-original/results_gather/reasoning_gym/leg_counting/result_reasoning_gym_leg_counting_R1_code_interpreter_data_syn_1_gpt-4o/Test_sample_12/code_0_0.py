# Define the number of legs for each type of animal
legs_per_animal = {
    'crabs': 10,
    'lobsters': 10,
    'humans': 2,
    'cows': 4,
    'bees': 6,
    'elephants': 4,
    'dogs': 4,
    'snakes': 0,
    'shrimps': 10
}

# Define the number of each type of animal
animal_counts = {
    'crabs': 2,
    'lobsters': 10,
    'humans': 1,
    'cows': 2,
    'bees': 3,
    'elephants': 13,
    'dogs': 9,
    'snakes': 12,
    'shrimps': 5
}

# Calculate the total number of legs
total_legs = sum(animal_counts[animal] * legs_per_animal[animal] for animal in animal_counts)

# Output the result
print(total_legs)