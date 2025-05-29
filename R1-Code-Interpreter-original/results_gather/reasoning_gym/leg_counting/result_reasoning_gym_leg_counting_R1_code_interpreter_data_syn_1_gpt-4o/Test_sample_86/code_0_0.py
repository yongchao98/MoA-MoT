# Define the number of legs for each type of animal
legs_per_animal = {
    "starfish": 0,
    "bird": 2,
    "beetle": 6,
    "sheep": 4,
    "sea_slug": 0
}

# Define the number of each type of animal
animal_counts = {
    "starfish": 5,
    "bird": 6,
    "beetle": 15,
    "sheep": 1,
    "sea_slug": 10
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

# Output the result
print(total_legs)