# Initial gifts
gifts = {
    "Alice": "orange ball",
    "Bob": "brown present",
    "Claire": "pink ball",
    "Dave": "blue present",
    "Eve": "green present",
    "Fred": "yellow present",
    "Gertrude": "white present"
}

# Swaps
swaps = [
    ("Dave", "Gertrude"),
    ("Gertrude", "Alice"),
    ("Claire", "Bob"),
    ("Eve", "Claire"),
    ("Fred", "Alice"),
    ("Gertrude", "Alice"),
    ("Bob", "Gertrude")
]

# Perform swaps
for person1, person2 in swaps:
    gifts[person1], gifts[person2] = gifts[person2], gifts[person1]

# Output the final gift of Gertrude
print(gifts["Gertrude"])