# Initial ball assignments
balls = {
    "Alice": "orange",
    "Bob": "green",
    "Claire": "blue",
    "Dave": "pink",
    "Eve": "purple"
}

# Sequence of swaps
swaps = [
    ("Dave", "Alice"),
    ("Bob", "Claire"),
    ("Bob", "Alice"),
    ("Alice", "Claire"),
    ("Eve", "Bob")
]

# Perform the swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Bob has at the end
print(balls["Bob"])