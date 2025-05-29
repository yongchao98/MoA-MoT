# Initial ball distribution
balls = {
    "Alice": "green",
    "Bob": "brown",
    "Claire": "purple",
    "Dave": "pink",
    "Eve": "white",
    "Fred": "red",
    "Gertrude": "orange"
}

# List of swaps
swaps = [
    ("Eve", "Dave"),
    ("Claire", "Eve"),
    ("Dave", "Claire"),
    ("Bob", "Fred"),
    ("Fred", "Dave"),
    ("Eve", "Alice"),
    ("Dave", "Gertrude")
]

# Perform the swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Alice has at the end
print(balls["Alice"])