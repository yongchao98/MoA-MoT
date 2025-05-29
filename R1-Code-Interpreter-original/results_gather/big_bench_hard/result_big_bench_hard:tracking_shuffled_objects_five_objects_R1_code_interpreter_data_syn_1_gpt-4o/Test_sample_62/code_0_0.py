# Initial ball positions
balls = {
    "Alice": "orange",
    "Bob": "black",
    "Claire": "green",
    "Dave": "brown",
    "Eve": "yellow"
}

# Swap operations
swaps = [
    ("Eve", "Alice"),
    ("Claire", "Bob"),
    ("Dave", "Bob"),
    ("Alice", "Eve"),
    ("Dave", "Claire")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the final ball Alice has
print(balls["Alice"])