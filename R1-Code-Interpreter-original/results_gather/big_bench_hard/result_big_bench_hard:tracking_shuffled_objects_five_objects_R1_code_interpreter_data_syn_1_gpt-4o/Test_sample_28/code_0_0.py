# Initial ball assignments
balls = {
    "Alice": "Purple",
    "Bob": "Black",
    "Claire": "Red",
    "Dave": "Yellow",
    "Eve": "Brown"
}

# Swaps
swaps = [
    ("Alice", "Dave"),
    ("Dave", "Eve"),
    ("Bob", "Dave"),
    ("Dave", "Claire"),
    ("Eve", "Claire")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Alice has at the end
print(balls["Alice"])