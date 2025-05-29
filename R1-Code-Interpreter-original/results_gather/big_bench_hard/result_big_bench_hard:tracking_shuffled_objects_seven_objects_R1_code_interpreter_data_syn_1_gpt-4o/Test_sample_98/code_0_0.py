# Initial ball distribution
balls = {
    "Alice": "green",
    "Bob": "white",
    "Claire": "yellow",
    "Dave": "pink",
    "Eve": "orange",
    "Fred": "black",
    "Gertrude": "brown"
}

# Swaps
swaps = [
    ("Bob", "Gertrude"),
    ("Fred", "Claire"),
    ("Dave", "Gertrude"),
    ("Bob", "Gertrude"),
    ("Alice", "Claire"),
    ("Gertrude", "Claire"),
    ("Eve", "Claire")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Fred has at the end
print(balls["Fred"])