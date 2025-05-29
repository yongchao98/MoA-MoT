# Initial pairings
pairings = {
    "Alice": "Karl",
    "Bob": "Jamie",
    "Claire": "Melissa",
    "Dave": "Lola",
    "Eve": "Ophelia"
}

# Function to switch partners
def switch_partners(dancer1, dancer2):
    partner1 = pairings[dancer1]
    partner2 = pairings[dancer2]
    pairings[dancer1] = partner2
    pairings[dancer2] = partner1

# Apply the partner switches
switch_partners("Bob", "Claire")
switch_partners("Alice", "Dave")
switch_partners("Bob", "Alice")
switch_partners("Claire", "Eve")
switch_partners("Eve", "Alice")

# Output the final partner of Claire
print(pairings["Claire"])