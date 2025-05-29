# Initial partners
partners = {
    "Alice": "Helga",
    "Bob": "Izzi",
    "Claire": "Ophelia",
    "Dave": "Melissa",
    "Eve": "Patrick",
    "Fred": "Rodrigo",
    "Gertrude": "Sam"
}

# Function to switch partners
def switch(partners, dancer1, dancer2):
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Perform the switches
switch(partners, "Alice", "Bob")
switch(partners, "Eve", "Gertrude")
switch(partners, "Alice", "Claire")
switch(partners, "Eve", "Dave")
switch(partners, "Bob", "Alice")
switch(partners, "Bob", "Fred")
switch(partners, "Bob", "Eve")

# Output Bob's final partner
print(partners["Bob"])