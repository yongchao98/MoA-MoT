# Initial pairings
partners = {
    'Alice': 'Izzi',
    'Bob': 'Jamie',
    'Claire': 'Ophelia',
    'Dave': 'Rodrigo',
    'Eve': 'Patrick',
    'Fred': 'Sam',
    'Gertrude': 'Lola'
}

# Function to switch partners
def switch(partners, dancer1, dancer2):
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Apply the switches
switch(partners, 'Eve', 'Claire')
switch(partners, 'Fred', 'Gertrude')
switch(partners, 'Eve', 'Bob')
switch(partners, 'Eve', 'Alice')
switch(partners, 'Dave', 'Bob')
switch(partners, 'Dave', 'Eve')
switch(partners, 'Alice', 'Dave')

# Find Gertrude's partner
gertrude_partner = partners['Gertrude']
print(gertrude_partner)