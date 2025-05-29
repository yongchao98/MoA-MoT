# Define the number of legs for each type of animal
legs_per_animal = {
    'insects': 6,
    'horses': 4,
    'giraffes': 4,
    'shrimps': 10,
    'humans': 2,
    'sheeps': 4,
    'lobsters': 10
}

# Define the number of each type of animal
animal_counts = {
    'insects': 9,
    'horses': 10,
    'giraffes': 4,
    'shrimps': 3,
    'humans': 9,
    'sheeps': 11,
    'lobsters': 11
}

# Calculate the total number of legs
total_legs = sum(legs_per_animal[animal] * animal_counts[animal] for animal in animal_counts)

# Output the result
print(total_legs)