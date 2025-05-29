# Number of legs per animal
legs_per_animal = {
    "sheep": 4,
    "dog": 4,
    "tiger": 4,
    "elephant": 4,
    "firefly": 6,
    "praying_mantis": 6,
    "cockroach": 6
}

# Number of each animal
animal_counts = {
    "sheep": 6,
    "dog": 15,
    "tiger": 1,
    "elephant": 13,
    "firefly": 15,
    "praying_mantis": 14,
    "cockroach": 2
}

# Calculate total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

print(total_legs)