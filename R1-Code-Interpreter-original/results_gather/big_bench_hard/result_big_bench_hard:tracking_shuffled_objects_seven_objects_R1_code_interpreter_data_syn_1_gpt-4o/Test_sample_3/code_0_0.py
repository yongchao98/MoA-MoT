# Initial pairings
partners = {
    "Alice": "Ophelia",
    "Bob": "Melissa",
    "Claire": "Jamie",
    "Dave": "Sam",
    "Eve": "Patrick",
    "Fred": "Rodrigo",
    "Gertrude": "Karl"
}

# Function to switch partners
def switch_partners(person1, person2):
    partners[person1], partners[person2] = partners[person2], partners[person1]

# Apply the partner switches
switch_partners("Dave", "Claire")
switch_partners("Alice", "Eve")
switch_partners("Eve", "Bob")
switch_partners("Claire", "Bob")
switch_partners("Fred", "Eve")
switch_partners("Gertrude", "Dave")
switch_partners("Dave", "Alice")

# Output Dave's final partner
print(partners["Dave"])