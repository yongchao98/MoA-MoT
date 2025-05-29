# Initial ball distribution
balls = {
    "Alice": "white",
    "Bob": "brown",
    "Claire": "blue",
    "Dave": "pink",
    "Eve": "purple",
    "Fred": "black",
    "Gertrude": "yellow"
}

# Swaps
swaps = [
    ("Fred", "Dave"),
    ("Gertrude", "Claire"),
    ("Dave", "Bob"),
    ("Claire", "Gertrude"),
    ("Alice", "Fred"),
    ("Bob", "Alice"),
    ("Eve", "Bob")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Gertrude has at the end
print(balls["Gertrude"])