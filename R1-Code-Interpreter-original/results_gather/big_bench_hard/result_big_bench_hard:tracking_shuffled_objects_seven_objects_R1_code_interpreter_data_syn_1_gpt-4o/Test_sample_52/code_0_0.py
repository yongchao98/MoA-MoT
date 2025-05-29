# Initial gifts
gifts = {
    "Alice": "blue present",
    "Bob": "orange ball",
    "Claire": "brown present",
    "Dave": "pink ball",
    "Eve": "black ball",
    "Fred": "green present",
    "Gertrude": "purple present"
}

# Swaps
swaps = [
    ("Fred", "Eve"),
    ("Bob", "Alice"),
    ("Alice", "Gertrude"),
    ("Dave", "Bob"),
    ("Dave", "Alice"),
    ("Bob", "Gertrude"),
    ("Claire", "Bob")
]

# Perform swaps
for person1, person2 in swaps:
    gifts[person1], gifts[person2] = gifts[person2], gifts[person1]

# Output the gift Dave has at the end
print(gifts["Dave"])