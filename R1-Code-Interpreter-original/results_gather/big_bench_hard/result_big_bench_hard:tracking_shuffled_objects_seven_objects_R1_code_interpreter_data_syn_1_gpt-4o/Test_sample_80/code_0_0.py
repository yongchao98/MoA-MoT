# Initial ball distribution
balls = {
    "Alice": "purple",
    "Bob": "blue",
    "Claire": "yellow",
    "Dave": "brown",
    "Eve": "black",
    "Fred": "white",
    "Gertrude": "orange"
}

# Swaps
swaps = [
    ("Dave", "Fred"),
    ("Claire", "Dave"),
    ("Eve", "Bob"),
    ("Gertrude", "Claire"),
    ("Alice", "Gertrude"),
    ("Bob", "Claire"),
    ("Bob", "Alice")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Eve has at the end
print(balls["Eve"])