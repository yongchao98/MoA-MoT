# Initial positions
positions = {
    "Alice": "left midfielder",
    "Bob": "right midfielder",
    "Claire": "striker",
    "Dave": "fullback",
    "Eve": "center midfielder",
    "Fred": "benchwarmer",
    "Gertrude": "right winger"
}

# Swaps
swaps = [
    ("Gertrude", "Bob"),
    ("Fred", "Claire"),
    ("Fred", "Dave"),
    ("Eve", "Dave"),
    ("Claire", "Alice"),
    ("Dave", "Bob"),
    ("Gertrude", "Alice")
]

# Perform swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Alice's final position
print(positions["Alice"])