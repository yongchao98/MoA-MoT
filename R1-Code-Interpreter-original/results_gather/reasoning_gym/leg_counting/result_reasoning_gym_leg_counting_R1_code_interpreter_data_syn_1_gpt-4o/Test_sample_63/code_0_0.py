# Define the number of legs for each type of animal
legs_per_animal = {
    'spider': 8,
    'sea_slug': 0,
    'duck': 2,
    'snake': 0,
    'chicken': 2
}

# Define the number of each type of animal
animal_counts = {
    'spider': 11,
    'sea_slug': 5,
    'duck': 13,
    'snake': 4,
    'chicken': 1
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

# Output the result
print(total_legs)