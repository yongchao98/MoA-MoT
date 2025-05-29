# Initial pairings
partners = {
    "Alice": "Jamie",
    "Bob": "Rodrigo",
    "Claire": "Patrick",
    "Dave": "Izzi",
    "Eve": "Melissa"
}

# Function to switch partners
def switch(partners, person1, person2):
    partners[person1], partners[person2] = partners[person2], partners[person1]

# Apply the partner switches
switch(partners, "Dave", "Eve")
switch(partners, "Bob", "Claire")
switch(partners, "Bob", "Eve")
switch(partners, "Claire", "Alice")
switch(partners, "Claire", "Bob")

# Output Dave's partner at the end
print(partners["Dave"])