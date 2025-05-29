# Initial positions
positions = {
    "Alice": "fullback",
    "Bob": "center midfielder",
    "Claire": "right midfielder",
    "Dave": "left winger",
    "Eve": "goalkeeper"
}

# Swaps
swaps = [
    ("Dave", "Eve"),
    ("Claire", "Bob"),
    ("Alice", "Bob"),
    ("Eve", "Dave"),
    ("Bob", "Alice")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Final position of Eve
print(positions["Eve"])