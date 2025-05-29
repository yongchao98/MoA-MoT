# Number of legs for each type of animal
legs_per_animal = {
    "horses": 4,
    "ants": 6,
    "cows": 4,
    "snakes": 0,
    "crabs": 10,
    "shrimps": 10,
    "spiders": 8,
    "elephants": 4,
    "leeches": 0
}

# Number of each type of animal
animal_counts = {
    "horses": 4,
    "ants": 2,
    "cows": 12,
    "snakes": 3,
    "crabs": 10,
    "shrimps": 3,
    "spiders": 8,
    "elephants": 7,
    "leeches": 6
}

# Calculate total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

print(total_legs)