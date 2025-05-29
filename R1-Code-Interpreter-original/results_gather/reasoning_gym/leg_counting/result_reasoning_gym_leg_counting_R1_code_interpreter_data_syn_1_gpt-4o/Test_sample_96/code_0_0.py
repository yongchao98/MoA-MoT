# Number of legs for each type of animal
legs_per_animal = {
    "wasp": 6,
    "scorpion": 8,
    "beetle": 6,
    "duck": 2,
    "firefly": 6,
    "sea_slug": 0
}

# Number of each type of animal
animal_counts = {
    "wasp": 3,
    "scorpion": 6,
    "beetle": 11,
    "duck": 9,
    "firefly": 5,
    "sea_slug": 11
}

# Calculate total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

print(total_legs)