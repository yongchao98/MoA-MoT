# Initial gifts
gifts = {
    "Alice": "brown present",
    "Bob": "purple present",
    "Claire": "orange ball",
    "Dave": "pink ball",
    "Eve": "blue present"
}

# Swap operations
swaps = [
    ("Claire", "Eve"),
    ("Claire", "Bob"),
    ("Alice", "Dave"),
    ("Bob", "Claire"),
    ("Bob", "Alice")
]

# Perform swaps
for person1, person2 in swaps:
    gifts[person1], gifts[person2] = gifts[person2], gifts[person1]

# Output the final gift with Eve
print(gifts["Eve"])