# Initial positions
positions = {
    "Alice": "striker",
    "Bob": "right winger",
    "Claire": "left winger",
    "Dave": "benchwarmer",
    "Eve": "goalkeeper",
    "Fred": "center midfielder",
    "Gertrude": "cheerleader"
}

# Swaps
swaps = [
    ("Eve", "Claire"),
    ("Gertrude", "Alice"),
    ("Fred", "Bob"),
    ("Dave", "Fred"),
    ("Fred", "Bob"),
    ("Bob", "Eve"),
    ("Claire", "Alice")
]

# Perform swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Gertrude's final position
print(positions["Gertrude"])