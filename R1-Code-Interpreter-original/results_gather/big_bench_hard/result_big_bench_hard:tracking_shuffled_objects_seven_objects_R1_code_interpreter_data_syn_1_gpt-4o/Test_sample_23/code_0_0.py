# Initial ball distribution
balls = {
    "Alice": "black",
    "Bob": "green",
    "Claire": "white",
    "Dave": "red",
    "Eve": "brown",
    "Fred": "purple",
    "Gertrude": "orange"
}

# Swaps
swaps = [
    ("Bob", "Claire"),
    ("Dave", "Eve"),
    ("Fred", "Bob"),
    ("Dave", "Gertrude"),
    ("Fred", "Bob"),
    ("Gertrude", "Fred"),
    ("Alice", "Fred")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Eve has at the end
print(balls["Eve"])