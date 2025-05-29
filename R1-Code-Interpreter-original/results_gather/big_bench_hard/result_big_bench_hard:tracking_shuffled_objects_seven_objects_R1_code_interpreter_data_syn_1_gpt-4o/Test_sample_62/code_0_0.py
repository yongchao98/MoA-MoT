# Initial pairings
pairings = {
    "Alice": "Melissa",
    "Bob": "Rodrigo",
    "Claire": "Karl",
    "Dave": "Ophelia",
    "Eve": "Lola",
    "Fred": "Jamie",
    "Gertrude": "Patrick"
}

# List of switches
switches = [
    ("Alice", "Eve"),
    ("Gertrude", "Claire"),
    ("Fred", "Claire"),
    ("Dave", "Alice"),
    ("Dave", "Fred"),
    ("Alice", "Bob"),
    ("Fred", "Bob")
]

# Apply each switch
for dancer1, dancer2 in switches:
    # Swap partners
    pairings[dancer1], pairings[dancer2] = pairings[dancer2], pairings[dancer1]

# Output the final partner of Gertrude
print(pairings["Gertrude"])