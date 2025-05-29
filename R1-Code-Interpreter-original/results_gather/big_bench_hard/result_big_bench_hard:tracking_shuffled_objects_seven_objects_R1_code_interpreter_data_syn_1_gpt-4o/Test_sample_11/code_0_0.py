# Initial pairings
partners = {
    "Alice": "Helga",
    "Bob": "Karl",
    "Claire": "Melissa",
    "Dave": "Izzi",
    "Eve": "Jamie",
    "Fred": "Patrick",
    "Gertrude": "Rodrigo"
}

# List of partner switches
switches = [
    ("Fred", "Claire"),
    ("Alice", "Bob"),
    ("Dave", "Eve"),
    ("Gertrude", "Fred"),
    ("Alice", "Bob"),
    ("Alice", "Fred"),
    ("Fred", "Claire")
]

# Apply each switch
for dancer1, dancer2 in switches:
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Output the final partner of Alice
print(partners["Alice"])