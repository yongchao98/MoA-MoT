# Initial ball assignments
balls = {
    "Alice": "orange",
    "Bob": "black",
    "Claire": "red",
    "Dave": "pink",
    "Eve": "blue"
}

# Swaps
swaps = [
    ("Eve", "Bob"),
    ("Claire", "Eve"),
    ("Bob", "Dave"),
    ("Eve", "Alice"),
    ("Dave", "Bob")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Alice has at the end
print(balls["Alice"])