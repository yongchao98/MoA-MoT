# Initial positions
positions = {
    "Alice": "left winger",
    "Bob": "goalkeeper",
    "Claire": "center midfielder",
    "Dave": "right midfielder",
    "Eve": "fullback",
    "Fred": "left midfielder",
    "Gertrude": "cheerleader"
}

# Swaps
swaps = [
    ("Gertrude", "Eve"),
    ("Gertrude", "Alice"),
    ("Eve", "Claire"),
    ("Fred", "Eve"),
    ("Claire", "Gertrude"),
    ("Claire", "Dave"),
    ("Fred", "Bob")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Gertrude's final position
print(positions["Gertrude"])