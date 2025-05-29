# Initial ball distribution
balls = {
    "Alice": "Purple",
    "Bob": "Green",
    "Claire": "Pink",
    "Dave": "Red",
    "Eve": "Blue",
    "Fred": "White",
    "Gertrude": "Black"
}

# Swaps
swaps = [
    ("Eve", "Dave"),
    ("Eve", "Fred"),
    ("Eve", "Dave"),
    ("Eve", "Alice"),
    ("Bob", "Alice"),
    ("Gertrude", "Claire"),
    ("Alice", "Dave")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the final ball Dave has
print(balls["Dave"])