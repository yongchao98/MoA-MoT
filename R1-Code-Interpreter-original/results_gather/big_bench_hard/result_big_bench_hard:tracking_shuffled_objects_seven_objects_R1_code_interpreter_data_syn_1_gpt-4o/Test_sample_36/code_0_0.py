# Initial positions
positions = {
    "Alice": "center midfielder",
    "Bob": "striker",
    "Claire": "left winger",
    "Dave": "cheerleader",
    "Eve": "fullback",
    "Fred": "right winger",
    "Gertrude": "benchwarmer"
}

# Swaps
swaps = [
    ("Alice", "Claire"),
    ("Gertrude", "Bob"),
    ("Eve", "Alice"),
    ("Claire", "Fred"),
    ("Gertrude", "Dave"),
    ("Eve", "Fred"),
    ("Fred", "Bob")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Dave's final position
print(positions["Dave"])