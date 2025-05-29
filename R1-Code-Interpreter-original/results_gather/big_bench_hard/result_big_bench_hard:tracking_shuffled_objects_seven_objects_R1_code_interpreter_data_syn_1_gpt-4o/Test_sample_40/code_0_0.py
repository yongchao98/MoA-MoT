# Initial pairings
partners = {
    "Alice": "Izzi",
    "Bob": "Melissa",
    "Claire": "Jamie",
    "Dave": "Rodrigo",
    "Eve": "Sam",
    "Fred": "Karl",
    "Gertrude": "Lola"
}

# List of switches
switches = [
    ("Dave", "Claire"),
    ("Claire", "Gertrude"),
    ("Dave", "Bob"),
    ("Fred", "Eve"),
    ("Dave", "Fred"),
    ("Claire", "Fred"),
    ("Alice", "Bob")
]

# Apply each switch
for dancer1, dancer2 in switches:
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Output the final partner of Alice
print(partners["Alice"])