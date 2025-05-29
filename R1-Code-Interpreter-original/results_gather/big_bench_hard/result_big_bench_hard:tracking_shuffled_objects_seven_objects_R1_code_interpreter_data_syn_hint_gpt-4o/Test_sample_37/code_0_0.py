# Initial ball distribution
balls = {
    "Alice": "pink",
    "Bob": "white",
    "Claire": "brown",
    "Dave": "blue",
    "Eve": "yellow",
    "Fred": "black",
    "Gertrude": "purple"
}

# Swaps
swaps = [
    ("Fred", "Gertrude"),
    ("Gertrude", "Bob"),
    ("Fred", "Claire"),
    ("Claire", "Dave"),
    ("Claire", "Eve"),
    ("Alice", "Bob"),
    ("Claire", "Eve")
]

# Perform swaps
for person1, person2 in swaps:
    balls[person1], balls[person2] = balls[person2], balls[person1]

# Output the ball Eve has at the end
print(balls["Eve"])