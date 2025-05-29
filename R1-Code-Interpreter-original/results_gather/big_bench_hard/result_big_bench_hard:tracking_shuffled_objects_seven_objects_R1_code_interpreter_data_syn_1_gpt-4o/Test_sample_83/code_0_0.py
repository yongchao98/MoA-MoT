# Initial pairings
pairings = {
    "Alice": "Ophelia",
    "Bob": "Melissa",
    "Claire": "Jamie",
    "Dave": "Sam",
    "Eve": "Patrick",
    "Fred": "Rodrigo",
    "Gertrude": "Karl"
}

# List of switches
switches = [
    ("Dave", "Claire"),
    ("Alice", "Eve"),
    ("Eve", "Bob"),
    ("Claire", "Bob"),
    ("Fred", "Eve"),
    ("Gertrude", "Dave"),
    ("Dave", "Alice")
]

# Apply each switch
for dancer1, dancer2 in switches:
    # Swap partners
    pairings[dancer1], pairings[dancer2] = pairings[dancer2], pairings[dancer1]

# Output Fred's partner at the end
print(pairings["Fred"])