# Initial ball distribution
balls = {
    "Alice": "pink",
    "Bob": "brown",
    "Claire": "black",
    "Dave": "orange",
    "Eve": "yellow"
}

# Swap operations
swaps = [
    ("Claire", "Alice"),
    ("Dave", "Bob"),
    ("Eve", "Bob"),
    ("Bob", "Alice"),
    ("Eve", "Claire")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Alice has at the end
print(balls["Alice"])