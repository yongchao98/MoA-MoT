# Initial ball distribution
balls = {
    "Alice": "Blue",
    "Bob": "Yellow",
    "Claire": "Pink",
    "Dave": "Black",
    "Eve": "White",
    "Fred": "Brown",
    "Gertrude": "Orange"
}

# Swaps
swaps = [
    ("Claire", "Gertrude"),
    ("Eve", "Alice"),
    ("Gertrude", "Claire"),
    ("Alice", "Dave"),
    ("Fred", "Dave"),
    ("Claire", "Dave"),
    ("Bob", "Eve")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Gertrude has at the end
print(balls["Gertrude"])