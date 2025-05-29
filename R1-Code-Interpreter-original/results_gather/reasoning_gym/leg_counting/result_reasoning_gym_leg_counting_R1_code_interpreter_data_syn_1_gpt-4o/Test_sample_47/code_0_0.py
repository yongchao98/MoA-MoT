# Define the number of legs for each type of animal
legs_per_animal = {
    "cockroach": 6,
    "butterfly": 6,
    "cat": 4,
    "giraffe": 4,
    "spider": 8
}

# Define the number of each type of animal
animal_counts = {
    "cockroach": 3,
    "butterfly": 2,
    "cat": 6,
    "giraffe": 1,
    "spider": 4
}

# Calculate the total number of legs
total_legs = sum(animal_counts[animal] * legs_per_animal[animal] for animal in animal_counts)

# Output the result
print(total_legs)