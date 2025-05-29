# Initial pairings
partners = {
    "Alice": "Melissa",
    "Bob": "Sam",
    "Claire": "Jamie",
    "Dave": "Izzi",
    "Eve": "Lola",
    "Fred": "Rodrigo",
    "Gertrude": "Ophelia"
}

# Partner switches
switches = [
    ("Gertrude", "Alice"),
    ("Dave", "Claire"),
    ("Bob", "Claire"),
    ("Gertrude", "Eve"),
    ("Claire", "Fred"),
    ("Fred", "Alice"),
    ("Gertrude", "Dave")
]

# Apply each switch
for dancer1, dancer2 in switches:
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Output Bob's partner at the end
print(partners["Bob"])