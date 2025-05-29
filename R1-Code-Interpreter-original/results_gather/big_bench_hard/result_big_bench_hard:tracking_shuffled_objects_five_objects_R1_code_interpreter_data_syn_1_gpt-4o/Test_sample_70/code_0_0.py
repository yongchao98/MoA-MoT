# Initial positions
positions = {
    "Alice": "cheerleader",
    "Bob": "goalkeeper",
    "Claire": "striker",
    "Dave": "left midfielder",
    "Eve": "right winger"
}

# Swaps
swaps = [
    ("Claire", "Eve"),
    ("Claire", "Dave"),
    ("Eve", "Claire"),
    ("Claire", "Alice"),
    ("Bob", "Dave")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Final position of Dave
print(positions["Dave"])