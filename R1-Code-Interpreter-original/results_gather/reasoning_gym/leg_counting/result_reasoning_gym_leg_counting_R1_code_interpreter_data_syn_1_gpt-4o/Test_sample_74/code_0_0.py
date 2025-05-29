# Define the number of legs for each animal
legs_per_animal = {
    "sea_slug": 0,
    "shrimp": 10,
    "crab": 10,
    "firefly": 6,
    "flatworm": 0,
    "tiger": 4,
    "wasp": 6,
    "jellyfish": 0,
    "bird": 2,
    "lobster": 10
}

# Define the number of each animal
animal_counts = {
    "sea_slug": 8,
    "shrimp": 6,
    "crab": 7,
    "firefly": 8,
    "flatworm": 12,
    "tiger": 11,
    "wasp": 5,
    "jellyfish": 12,
    "bird": 4,
    "lobster": 12
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

# Output the total number of legs
print(total_legs)