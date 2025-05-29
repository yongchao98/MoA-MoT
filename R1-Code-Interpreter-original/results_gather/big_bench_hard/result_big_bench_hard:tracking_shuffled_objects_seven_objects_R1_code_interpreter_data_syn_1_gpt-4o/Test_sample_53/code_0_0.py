# Initial ball distribution
balls = {
    "Alice": "green",
    "Bob": "orange",
    "Claire": "white",
    "Dave": "pink",
    "Eve": "blue",
    "Fred": "yellow",
    "Gertrude": "black"
}

# Swaps
swaps = [
    ("Eve", "Gertrude"),
    ("Fred", "Bob"),
    ("Fred", "Claire"),
    ("Eve", "Fred"),
    ("Fred", "Claire"),
    ("Alice", "Bob"),
    ("Dave", "Gertrude")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Gertrude has at the end
print(balls["Gertrude"])