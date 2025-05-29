# Initial pairings
partners = {
    'Alice': 'Lola',
    'Bob': 'Helga',
    'Claire': 'Sam',
    'Dave': 'Melissa',
    'Eve': 'Jamie',
    'Fred': 'Rodrigo',
    'Gertrude': 'Patrick'
}

# Function to switch partners
def switch(partners, person1, person2):
    partners[person1], partners[person2] = partners[person2], partners[person1]

# Apply the switches
switch(partners, 'Dave', 'Claire')
switch(partners, 'Gertrude', 'Claire')
switch(partners, 'Claire', 'Eve')
switch(partners, 'Gertrude', 'Eve')
switch(partners, 'Fred', 'Claire')
switch(partners, 'Bob', 'Alice')
switch(partners, 'Eve', 'Claire')

# Find Gertrude's partner
gertrude_partner = partners['Gertrude']
print(gertrude_partner)