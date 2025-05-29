# Initial gifts
gifts = {
    "Alice": "green present",
    "Bob": "orange ball",
    "Claire": "black ball",
    "Dave": "blue present",
    "Eve": "red present",
    "Fred": "brown present",
    "Gertrude": "pink ball"
}

# Swaps
swaps = [
    ("Fred", "Eve"),
    ("Bob", "Gertrude"),
    ("Alice", "Claire"),
    ("Dave", "Fred"),
    ("Alice", "Fred"),
    ("Dave", "Eve"),
    ("Alice", "Eve")
]

# Perform swaps
for person1, person2 in swaps:
    gifts[person1], gifts[person2] = gifts[person2], gifts[person1]

# Output the final gift Eve has
print(gifts["Eve"])