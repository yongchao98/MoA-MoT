# Initial gifts
gifts = {
    "Alice": "Blue",
    "Bob": "Brown",
    "Claire": "Orange",
    "Dave": "Red",
    "Eve": "Black",
    "Fred": "Purple",
    "Gertrude": "Green"
}

# Swaps
swaps = [
    ("Claire", "Gertrude"),
    ("Fred", "Gertrude"),
    ("Gertrude", "Eve"),
    ("Dave", "Eve"),
    ("Dave", "Alice"),
    ("Eve", "Alice"),
    ("Bob", "Dave")
]

# Perform swaps
for person1, person2 in swaps:
    gifts[person1], gifts[person2] = gifts[person2], gifts[person1]

# Output the final gift Dave has
print(gifts["Dave"])