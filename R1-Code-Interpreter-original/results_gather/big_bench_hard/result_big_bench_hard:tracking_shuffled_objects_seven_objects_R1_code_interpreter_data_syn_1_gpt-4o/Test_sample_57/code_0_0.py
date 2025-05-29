# Initial gifts
gifts = {
    "Alice": "green present",
    "Bob": "white present",
    "Claire": "orange ball",
    "Dave": "blue present",
    "Eve": "red present",
    "Fred": "black ball",
    "Gertrude": "pink ball"
}

# Swaps
swaps = [
    ("Eve", "Dave"),
    ("Fred", "Gertrude"),
    ("Fred", "Eve"),
    ("Claire", "Eve"),
    ("Claire", "Fred"),
    ("Bob", "Alice"),
    ("Bob", "Claire")
]

# Perform swaps
for person1, person2 in swaps:
    gifts[person1], gifts[person2] = gifts[person2], gifts[person1]

# Output the gift Alice has at the end
print(gifts["Alice"])