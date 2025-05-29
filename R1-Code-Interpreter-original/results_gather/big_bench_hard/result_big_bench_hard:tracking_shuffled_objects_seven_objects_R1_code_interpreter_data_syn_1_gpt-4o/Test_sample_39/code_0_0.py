# Initial pairings
partners = {
    'Alice': 'Jamie',
    'Bob': 'Ophelia',
    'Claire': 'Lola',
    'Dave': 'Patrick',
    'Eve': 'Rodrigo',
    'Fred': 'Melissa',
    'Gertrude': 'Karl'
}

# Function to switch partners
def switch(partners, dancer1, dancer2):
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Apply the switches
switch(partners, 'Eve', 'Dave')
switch(partners, 'Fred', 'Eve')
switch(partners, 'Alice', 'Gertrude')
switch(partners, 'Dave', 'Bob')
switch(partners, 'Alice', 'Claire')
switch(partners, 'Alice', 'Gertrude')
switch(partners, 'Bob', 'Claire')

# Find Gertrude's partner at the end
gertrude_partner = partners['Gertrude']
print(gertrude_partner)