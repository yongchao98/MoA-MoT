# Initial gifts
gifts = {
    "Alice": "pink ball",
    "Bob": "brown present",
    "Claire": "green present",
    "Dave": "red present",
    "Eve": "purple present"
}

# Swaps
swaps = [
    ("Claire", "Bob"),
    ("Claire", "Dave"),
    ("Bob", "Alice"),
    ("Alice", "Eve"),
    ("Dave", "Alice")
]

# Perform swaps
for person1, person2 in swaps:
    gifts[person1], gifts[person2] = gifts[person2], gifts[person1]

# Output the final gift Claire has
print(gifts["Claire"])