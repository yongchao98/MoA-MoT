# Initial pairings
partners = {
    "Alice": "Helga",
    "Bob": "Rodrigo",
    "Claire": "Melissa",
    "Dave": "Patrick",
    "Eve": "Jamie"
}

# Partner switches
switches = [
    ("Eve", "Bob"),
    ("Alice", "Bob"),
    ("Dave", "Eve"),
    ("Eve", "Claire"),
    ("Alice", "Claire")
]

# Function to switch partners
def switch_partners(partners, person1, person2):
    partners[person1], partners[person2] = partners[person2], partners[person1]

# Apply each switch
for person1, person2 in switches:
    switch_partners(partners, person1, person2)

# Output the final partner of Eve
print(partners["Eve"])