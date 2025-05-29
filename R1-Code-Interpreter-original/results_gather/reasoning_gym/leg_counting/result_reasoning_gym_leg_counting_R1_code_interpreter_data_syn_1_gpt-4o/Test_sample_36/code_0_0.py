# Number of legs for each animal
legs_per_animal = {
    "scorpions": 8,
    "humans": 2,
    "giraffes": 4,
    "crickets": 6,
    "fireflies": 6,
    "tigers": 4,
    "beetles": 6,
    "crabs": 10,
    "ants": 6,
    "cats": 4
}

# Number of each animal
animal_counts = {
    "scorpions": 8,
    "humans": 9,
    "giraffes": 6,
    "crickets": 3,
    "fireflies": 9,
    "tigers": 11,
    "beetles": 13,
    "crabs": 13,
    "ants": 13,
    "cats": 2
}

# Calculate total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

print(total_legs)