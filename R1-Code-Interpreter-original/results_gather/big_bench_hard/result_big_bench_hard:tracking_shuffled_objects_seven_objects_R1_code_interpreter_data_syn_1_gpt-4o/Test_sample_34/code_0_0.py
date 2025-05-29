# Initial positions
positions = {
    "Alice": "right midfielder",
    "Bob": "fullback",
    "Claire": "center midfielder",
    "Dave": "left midfielder",
    "Eve": "striker",
    "Fred": "cheerleader",
    "Gertrude": "goalkeeper"
}

# Swaps
swaps = [
    ("Gertrude", "Bob"),
    ("Eve", "Alice"),
    ("Gertrude", "Dave"),
    ("Dave", "Fred"),
    ("Claire", "Eve"),
    ("Fred", "Gertrude"),
    ("Fred", "Claire")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Bob's final position
print(positions["Bob"])