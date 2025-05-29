# Initial gifts
gifts = {
    "Alice": "orange",
    "Bob": "brown",
    "Claire": "pink",
    "Dave": "blue",
    "Eve": "green",
    "Fred": "yellow",
    "Gertrude": "white"
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

# Output the final gift that Bob has
print(gifts["Bob"])