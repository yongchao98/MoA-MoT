# Initial gifts
gifts = {
    "Alice": "Red",
    "Bob": "Purple",
    "Claire": "Green",
    "Dave": "Pink",
    "Eve": "White",
    "Fred": "Yellow",
    "Gertrude": "Brown"
}

# Swaps
swaps = [
    ("Dave", "Fred"),
    ("Gertrude", "Dave"),
    ("Dave", "Fred"),
    ("Claire", "Alice"),
    ("Claire", "Dave"),
    ("Bob", "Dave"),
    ("Alice", "Eve")
]

# Perform swaps
for person1, person2 in swaps:
    gifts[person1], gifts[person2] = gifts[person2], gifts[person1]

# Output the final gift Fred has
print(gifts["Fred"])