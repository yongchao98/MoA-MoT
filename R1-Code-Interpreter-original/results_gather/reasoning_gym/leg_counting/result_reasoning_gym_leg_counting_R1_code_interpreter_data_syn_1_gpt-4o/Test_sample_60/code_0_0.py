# Define the number of animals and their respective leg counts
animal_legs = {
    "crabs": (10, 10),
    "shrimps": (11, 10),
    "deers": (11, 4),
    "birds": (6, 2),
    "firefly": (1, 6),
    "cat": (1, 4),
    "grasshoppers": (14, 6),
    "scorpion": (1, 8),
    "jellyfish": (11, 0)
}

# Calculate the total number of legs
total_legs = sum(count * legs for count, legs in animal_legs.values())

# Output the result
print(total_legs)