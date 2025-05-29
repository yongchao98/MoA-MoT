# Initial ball assignments
balls = {
    "Alice": "blue",
    "Bob": "red",
    "Claire": "brown",
    "Dave": "purple",
    "Eve": "black",
    "Fred": "pink",
    "Gertrude": "green"
}

# List of swaps
swaps = [
    ("Fred", "Alice"),
    ("Claire", "Gertrude"),
    ("Bob", "Dave"),
    ("Eve", "Alice"),
    ("Claire", "Bob"),
    ("Dave", "Eve"),
    ("Claire", "Bob")
]

# Perform the swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Alice has at the end
print(balls["Alice"])