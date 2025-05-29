# Initial ball distribution
balls = {
    "Alice": "Blue",
    "Bob": "Green",
    "Claire": "Brown",
    "Dave": "Pink",
    "Eve": "Orange",
    "Fred": "Black",
    "Gertrude": "White"
}

# Sequence of swaps
swaps = [
    ("Gertrude", "Alice"),
    ("Claire", "Eve"),
    ("Bob", "Claire"),
    ("Claire", "Eve"),
    ("Eve", "Fred"),
    ("Bob", "Fred"),
    ("Dave", "Eve")
]

# Perform the swaps
for person1, person2 in swaps:
    balls[person1], balls[person2] = balls[person2], balls[person1]

# Output the ball Alice has at the end
print(balls["Alice"])