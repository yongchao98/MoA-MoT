# Initial presents
presents = {
    "Alice": "Green",
    "Bob": "Red",
    "Claire": "White",
    "Dave": "Blue",
    "Eve": "Yellow",
    "Fred": "Black",
    "Gertrude": "Purple"
}

# Swaps
swaps = [
    ("Bob", "Alice"),
    ("Dave", "Alice"),
    ("Gertrude", "Claire"),
    ("Dave", "Claire"),
    ("Alice", "Eve"),
    ("Gertrude", "Claire"),
    ("Fred", "Bob")
]

# Perform swaps
for person1, person2 in swaps:
    presents[person1], presents[person2] = presents[person2], presents[person1]

# Output the final present Alice has
print(presents["Alice"])