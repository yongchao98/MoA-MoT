# Define the number of legs for each type of animal
legs_per_animal = {
    "grasshopper": 6,
    "spider": 8,
    "tiger": 4,
    "chicken": 2,
    "starfish": 0,
    "ant": 6,
    "snake": 0
}

# Define the count of each type of animal
animal_counts = {
    "grasshopper": 2,
    "spider": 8,
    "tiger": 1,
    "chicken": 2,
    "starfish": 5,
    "ant": 13,
    "snake": 2
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

# Output the result
print(total_legs)