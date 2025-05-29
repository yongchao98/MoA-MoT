# Number of legs for each type of animal
legs_per_animal = {
    'ant': 6,
    'elephant': 4,
    'cat': 4,
    'bee': 6,
    'sheep': 4,
    'deer': 4,
    'horse': 4,
    'cricket': 6,
    'beetle': 6
}

# Number of each type of animal
animal_counts = {
    'ant': 13,
    'elephant': 10,
    'cat': 1,
    'bee': 13,
    'sheep': 4,
    'deer': 11,
    'horse': 14,
    'cricket': 1,
    'beetle': 8
}

# Calculate total number of legs
total_legs = sum(legs_per_animal[animal] * count for animal, count in animal_counts.items())

print(total_legs)