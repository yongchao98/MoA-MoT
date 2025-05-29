# Number of legs for each type of animal
legs_per_animal = {
    "chicken": 2,
    "firefly": 6,
    "starfish": 0,
    "butterfly": 6,
    "cockroach": 6,
    "lobster": 10
}

# Number of each type of animal
animal_counts = {
    "chicken": 5,
    "firefly": 8,
    "starfish": 5,
    "butterfly": 6,
    "cockroach": 2,
    "lobster": 9
}

# Calculate total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

print(total_legs)