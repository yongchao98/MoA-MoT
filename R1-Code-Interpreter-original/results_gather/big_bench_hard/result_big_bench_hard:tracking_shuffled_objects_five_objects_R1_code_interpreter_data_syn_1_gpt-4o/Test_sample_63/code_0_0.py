# Initial pairings
partners = {
    "Alice": "Sam",
    "Bob": "Karl",
    "Claire": "Izzi",
    "Dave": "Patrick",
    "Eve": "Helga"
}

# Partner switches
switches = [
    ("Bob", "Claire"),
    ("Eve", "Dave"),
    ("Bob", "Dave"),
    ("Alice", "Bob"),
    ("Bob", "Eve")
]

# Function to switch partners
def switch_partners(partners, person1, person2):
    partners[person1], partners[person2] = partners[person2], partners[person1]

# Apply each switch
for person1, person2 in switches:
    switch_partners(partners, person1, person2)

# Output who Bob is dancing with at the end
print(partners["Bob"])