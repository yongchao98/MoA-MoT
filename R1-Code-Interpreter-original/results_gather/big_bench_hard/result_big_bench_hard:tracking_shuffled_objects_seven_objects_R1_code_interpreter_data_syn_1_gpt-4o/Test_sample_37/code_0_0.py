# Initial ball distribution
balls = {
    "Alice": "Pink",
    "Bob": "White",
    "Claire": "Brown",
    "Dave": "Blue",
    "Eve": "Yellow",
    "Fred": "Black",
    "Gertrude": "Purple"
}

# List of swaps
swaps = [
    ("Fred", "Gertrude"),
    ("Gertrude", "Bob"),
    ("Fred", "Claire"),
    ("Claire", "Dave"),
    ("Claire", "Eve"),
    ("Alice", "Bob"),
    ("Claire", "Eve")
]

# Perform the swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Eve has at the end
print(balls["Eve"])