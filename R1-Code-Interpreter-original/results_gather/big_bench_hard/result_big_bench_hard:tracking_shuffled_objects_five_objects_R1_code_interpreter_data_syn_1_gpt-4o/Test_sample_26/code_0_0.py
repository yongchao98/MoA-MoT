# Initial gifts
gifts = {
    "Alice": "pink ball",
    "Bob": "red present",
    "Claire": "green present",
    "Dave": "orange ball",
    "Eve": "black ball"
}

# Swaps
swaps = [
    ("Alice", "Dave"),
    ("Dave", "Eve"),
    ("Dave", "Alice"),
    ("Bob", "Claire"),
    ("Bob", "Eve")
]

# Perform swaps
for person1, person2 in swaps:
    gifts[person1], gifts[person2] = gifts[person2], gifts[person1]

# Output the gift Bob has at the end
print(gifts["Bob"])