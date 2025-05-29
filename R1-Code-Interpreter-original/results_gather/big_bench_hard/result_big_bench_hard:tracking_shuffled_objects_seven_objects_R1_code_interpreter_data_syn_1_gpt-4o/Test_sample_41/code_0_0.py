# Initial presents
presents = {
    "Alice": "Red",
    "Bob": "Purple",
    "Claire": "White",
    "Dave": "Brown",
    "Eve": "Green",
    "Fred": "Blue",
    "Gertrude": "Yellow"
}

# Swaps
swaps = [
    ("Dave", "Bob"),
    ("Alice", "Bob"),
    ("Gertrude", "Alice"),
    ("Claire", "Alice"),
    ("Bob", "Dave"),
    ("Fred", "Alice"),
    ("Fred", "Eve")
]

# Perform swaps
for person1, person2 in swaps:
    presents[person1], presents[person2] = presents[person2], presents[person1]

# Output the present Gertrude has
print(presents["Gertrude"])