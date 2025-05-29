# Initial partners
partners = {
    "Alice": "Helga",
    "Bob": "Karl",
    "Claire": "Melissa",
    "Dave": "Ophelia",
    "Eve": "Sam"
}

# Partner switches
switches = [
    ("Alice", "Dave"),
    ("Eve", "Alice"),
    ("Bob", "Claire"),
    ("Alice", "Claire"),
    ("Dave", "Eve")
]

# Apply each switch
for dancer1, dancer2 in switches:
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Output the final partner of Dave
print(partners["Dave"])