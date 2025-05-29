# Initial positions
positions = {
    "Alice": "striker",
    "Bob": "goalkeeper",
    "Claire": "left winger",
    "Dave": "benchwarmer",
    "Eve": "center midfielder",
    "Fred": "fullback",
    "Gertrude": "cheerleader"
}

# Swaps
swaps = [
    ("Alice", "Fred"),
    ("Claire", "Gertrude"),
    ("Dave", "Alice"),
    ("Bob", "Dave"),
    ("Bob", "Eve"),
    ("Eve", "Alice"),
    ("Gertrude", "Dave")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Final position of Dave
print(positions["Dave"])