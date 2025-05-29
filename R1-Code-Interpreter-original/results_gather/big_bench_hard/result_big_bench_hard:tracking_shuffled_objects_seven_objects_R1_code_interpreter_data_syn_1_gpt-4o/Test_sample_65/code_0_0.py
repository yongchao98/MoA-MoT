# Initial ball distribution
balls = {
    "Alice": "Black",
    "Bob": "Blue",
    "Claire": "Red",
    "Dave": "Purple",
    "Eve": "Pink",
    "Fred": "Brown",
    "Gertrude": "White"
}

# Sequence of swaps
swaps = [
    ("Fred", "Gertrude"),
    ("Fred", "Alice"),
    ("Alice", "Claire"),
    ("Bob", "Dave"),
    ("Fred", "Eve"),
    ("Gertrude", "Fred"),
    ("Gertrude", "Bob")
]

# Perform the swaps
for person1, person2 in swaps:
    balls[person1], balls[person2] = balls[person2], balls[person1]

# Output the final ball Alice has
print(balls["Alice"])