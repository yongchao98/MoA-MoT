# Initial partners
partners = {
    "Alice": "Melissa",
    "Bob": "Sam",
    "Claire": "Izzi",
    "Dave": "Rodrigo",
    "Eve": "Karl",
    "Fred": "Lola",
    "Gertrude": "Patrick"
}

# Function to switch partners
def switch_partners(person1, person2):
    partners[person1], partners[person2] = partners[person2], partners[person1]

# Apply the partner switches
switch_partners("Dave", "Bob")
switch_partners("Claire", "Bob")
switch_partners("Eve", "Bob")
switch_partners("Gertrude", "Bob")
switch_partners("Claire", "Dave")
switch_partners("Alice", "Gertrude")
switch_partners("Bob", "Fred")

# Output Bob's final partner
print(partners["Bob"])