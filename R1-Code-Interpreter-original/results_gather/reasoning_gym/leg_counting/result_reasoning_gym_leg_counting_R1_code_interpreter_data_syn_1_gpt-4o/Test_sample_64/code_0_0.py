# Define the number of legs for each type of animal
legs_per_animal = {
    "insect": 6,
    "wasp": 6,
    "lion": 4,
    "bird": 2,
    "ant": 6,
    "jellyfish": 0,
    "butterfly": 6,
    "cockroach": 6
}

# Define the number of each type of animal
animal_counts = {
    "insect": 1,
    "wasp": 3,
    "lion": 7,
    "bird": 11,
    "ant": 6,
    "jellyfish": 2,
    "butterfly": 5,
    "cockroach": 6
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

# Output the total number of legs
print(total_legs)