# Initial partners
partners = {
    "Alice": "Helga",
    "Bob": "Jamie",
    "Claire": "Ophelia",
    "Dave": "Karl",
    "Eve": "Rodrigo"
}

# Partner switches
switches = [
    ("Eve", "Bob"),
    ("Dave", "Claire"),
    ("Claire", "Bob"),
    ("Dave", "Alice"),
    ("Claire", "Eve")
]

# Apply each switch
for dancer1, dancer2 in switches:
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Find Claire's partner at the end
claire_partner = partners["Claire"]
print(claire_partner)