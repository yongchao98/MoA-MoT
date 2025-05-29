# Initial ball distribution
balls = {
    "Alice": "orange",
    "Bob": "black",
    "Claire": "pink",
    "Dave": "white",
    "Eve": "red"
}

# Swaps
swaps = [
    ("Alice", "Dave"),
    ("Bob", "Claire"),
    ("Claire", "Dave"),
    ("Dave", "Bob"),
    ("Alice", "Eve")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Eve has at the end
print(balls["Eve"])