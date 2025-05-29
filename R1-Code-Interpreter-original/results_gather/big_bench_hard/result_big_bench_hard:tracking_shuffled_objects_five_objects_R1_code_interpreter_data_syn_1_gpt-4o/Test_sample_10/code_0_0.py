# Initial ball assignments
balls = {
    "Alice": "Red",
    "Bob": "Black",
    "Claire": "Blue",
    "Dave": "Yellow",
    "Eve": "Brown"
}

# Swap sequence
swaps = [
    ("Claire", "Eve"),
    ("Alice", "Bob"),
    ("Eve", "Dave"),
    ("Claire", "Dave"),
    ("Claire", "Alice")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the final ball Claire has
print(balls["Claire"])