# Initial positions
positions = {
    "Alice": "right winger",
    "Bob": "goalkeeper",
    "Claire": "benchwarmer",
    "Dave": "center midfielder",
    "Eve": "left midfielder"
}

# Swaps
swaps = [
    ("Claire", "Eve"),
    ("Claire", "Alice"),
    ("Alice", "Dave"),
    ("Bob", "Eve"),
    ("Alice", "Claire")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Claire's final position
print(positions["Claire"])