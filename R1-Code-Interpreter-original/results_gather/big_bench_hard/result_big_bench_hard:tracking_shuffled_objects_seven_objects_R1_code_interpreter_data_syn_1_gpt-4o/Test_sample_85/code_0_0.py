# Initial ball distribution
balls = {
    "Alice": "Blue",
    "Bob": "Purple",
    "Claire": "Orange",
    "Dave": "Pink",
    "Eve": "Red",
    "Fred": "Green",
    "Gertrude": "Brown"
}

# Swap operations
swaps = [
    ("Gertrude", "Eve"),
    ("Dave", "Fred"),
    ("Claire", "Fred"),
    ("Eve", "Bob"),
    ("Claire", "Alice"),
    ("Gertrude", "Eve"),
    ("Claire", "Alice")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Eve has at the end
print(balls["Eve"])