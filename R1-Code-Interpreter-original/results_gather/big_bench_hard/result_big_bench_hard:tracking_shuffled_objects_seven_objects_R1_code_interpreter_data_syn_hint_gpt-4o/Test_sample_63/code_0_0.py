# Initial ball positions
balls = {
    "Alice": "orange",
    "Bob": "pink",
    "Claire": "brown",
    "Dave": "green",
    "Eve": "yellow",
    "Fred": "red",
    "Gertrude": "blue"
}

# Swaps
swaps = [
    ("Bob", "Gertrude"),
    ("Dave", "Bob"),
    ("Eve", "Gertrude"),
    ("Bob", "Alice"),
    ("Eve", "Fred"),
    ("Claire", "Dave"),
    ("Gertrude", "Fred")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Bob has at the end
print(balls["Bob"])