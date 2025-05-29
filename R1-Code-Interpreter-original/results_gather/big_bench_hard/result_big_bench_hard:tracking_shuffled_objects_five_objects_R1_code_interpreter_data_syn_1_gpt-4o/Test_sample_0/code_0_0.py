# Initial pairings
partners = {
    "Alice": "Patrick",
    "Bob": "Sam",
    "Claire": "Jamie",
    "Dave": "Lola",
    "Eve": "Melissa"
}

# Function to switch partners
def switch_partners(person1, person2):
    partners[person1], partners[person2] = partners[person2], partners[person1]

# Apply the partner switches
switch_partners("Dave", "Eve")
switch_partners("Dave", "Alice")
switch_partners("Eve", "Alice")
switch_partners("Claire", "Bob")
switch_partners("Dave", "Alice")

# Output the final partner of Alice
print(partners["Alice"])