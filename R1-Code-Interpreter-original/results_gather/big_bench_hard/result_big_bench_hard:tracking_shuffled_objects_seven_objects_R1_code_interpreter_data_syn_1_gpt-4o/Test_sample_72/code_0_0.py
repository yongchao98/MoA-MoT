# Initial positions
positions = {
    "Alice": "left midfielder",
    "Bob": "center midfielder",
    "Claire": "benchwarmer",
    "Dave": "cheerleader",
    "Eve": "striker",
    "Fred": "right midfielder",
    "Gertrude": "right winger"
}

# Swaps
swaps = [
    ("Dave", "Gertrude"),
    ("Dave", "Fred"),
    ("Alice", "Dave"),
    ("Alice", "Gertrude"),
    ("Bob", "Eve"),
    ("Dave", "Alice"),
    ("Claire", "Dave")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Gertrude's final position
print(positions["Gertrude"])