# Initial positions
positions = {
    "Alice": "benchwarmer",
    "Bob": "striker",
    "Claire": "right winger",
    "Dave": "center midfielder",
    "Eve": "cheerleader",
    "Fred": "left winger",
    "Gertrude": "goalkeeper"
}

# List of swaps
swaps = [
    ("Eve", "Alice"),
    ("Alice", "Fred"),
    ("Alice", "Gertrude"),
    ("Bob", "Eve"),
    ("Alice", "Bob"),
    ("Eve", "Claire"),
    ("Alice", "Dave")
]

# Process each swap
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output the final position of Gertrude
print(positions["Gertrude"])