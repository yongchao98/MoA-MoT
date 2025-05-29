# Initial pairings
partners = {
    "Alice": "Sam",
    "Bob": "Helga",
    "Claire": "Ophelia",
    "Dave": "Lola",
    "Eve": "Melissa",
    "Fred": "Patrick",
    "Gertrude": "Jamie"
}

# List of switches
switches = [
    ("Alice", "Dave"),
    ("Claire", "Alice"),
    ("Fred", "Dave"),
    ("Eve", "Fred"),
    ("Claire", "Fred"),
    ("Alice", "Eve"),
    ("Gertrude", "Bob")
]

# Apply each switch
for dancer1, dancer2 in switches:
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Output the partner of Dave at the end
print(partners["Dave"])