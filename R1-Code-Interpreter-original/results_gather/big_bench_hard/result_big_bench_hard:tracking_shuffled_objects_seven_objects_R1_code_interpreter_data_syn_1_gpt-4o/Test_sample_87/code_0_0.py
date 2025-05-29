# Initial ball distribution
balls = {
    "Alice": "Blue",
    "Bob": "Pink",
    "Claire": "Brown",
    "Dave": "Yellow",
    "Eve": "Black",
    "Fred": "Red",
    "Gertrude": "Purple"
}

# List of swaps
swaps = [
    ("Claire", "Dave"),
    ("Eve", "Alice"),
    ("Gertrude", "Eve"),
    ("Claire", "Alice"),
    ("Dave", "Eve"),
    ("Alice", "Gertrude"),
    ("Fred", "Bob")
]

# Perform the swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Eve has at the end
print(balls["Eve"])