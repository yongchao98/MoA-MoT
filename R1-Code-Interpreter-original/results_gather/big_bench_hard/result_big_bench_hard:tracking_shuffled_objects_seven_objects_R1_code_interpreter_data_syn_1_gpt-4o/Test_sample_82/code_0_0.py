# Initial positions
positions = {
    "Alice": "fullback",
    "Bob": "benchwarmer",
    "Claire": "cheerleader",
    "Dave": "left midfielder",
    "Eve": "right winger",
    "Fred": "left winger",
    "Gertrude": "goalkeeper"
}

# Swaps
swaps = [
    ("Eve", "Bob"),
    ("Gertrude", "Bob"),
    ("Bob", "Dave"),
    ("Gertrude", "Dave"),
    ("Bob", "Claire"),
    ("Fred", "Bob"),
    ("Alice", "Eve")
]

# Perform swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Gertrude's final position
print(positions["Gertrude"])