# Initial positions
positions = {
    "Alice": "left midfielder",
    "Bob": "benchwarmer",
    "Claire": "right midfielder",
    "Dave": "cheerleader",
    "Eve": "left winger",
    "Fred": "goalkeeper",
    "Gertrude": "right winger"
}

# Swaps
swaps = [
    ("Fred", "Claire"),
    ("Fred", "Dave"),
    ("Alice", "Bob"),
    ("Gertrude", "Eve"),
    ("Dave", "Bob"),
    ("Dave", "Eve"),
    ("Dave", "Alice")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Gertrude's final position
print(positions["Gertrude"])