# Initial assignment
presents = {
    "Alice": "blue present",
    "Bob": "green present",
    "Claire": "brown present",
    "Dave": "purple present",
    "Eve": "orange ball",
    "Fred": "pink ball",
    "Gertrude": "red present"
}

# Swaps
swaps = [
    ("Eve", "Claire"),
    ("Fred", "Claire"),
    ("Dave", "Bob"),
    ("Alice", "Claire"),
    ("Dave", "Alice"),
    ("Gertrude", "Eve"),
    ("Claire", "Gertrude")
]

# Perform swaps
for person1, person2 in swaps:
    presents[person1], presents[person2] = presents[person2], presents[person1]

# Output Bob's present
print(presents["Bob"])