# Initial gifts
gifts = {
    "Alice": "black",
    "Bob": "orange",
    "Claire": "white",
    "Dave": "green",
    "Eve": "blue",
    "Fred": "red",
    "Gertrude": "yellow"
}

# List of swaps
swaps = [
    ("Alice", "Eve"),
    ("Alice", "Fred"),
    ("Gertrude", "Dave"),
    ("Fred", "Eve"),
    ("Claire", "Bob"),
    ("Claire", "Dave"),
    ("Dave", "Alice")
]

# Perform the swaps
for person1, person2 in swaps:
    gifts[person1], gifts[person2] = gifts[person2], gifts[person1]

# Output the final gift Alice has
print(gifts["Alice"])