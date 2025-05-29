# Initial pairings
partners = {
    'Alice': 'Ophelia',
    'Bob': 'Patrick',
    'Claire': 'Karl',
    'Dave': 'Rodrigo',
    'Eve': 'Izzi',
    'Fred': 'Helga',
    'Gertrude': 'Jamie'
}

# Function to switch partners
def switch(partners, dancer1, dancer2):
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Apply the switches
switch(partners, 'Claire', 'Eve')
switch(partners, 'Gertrude', 'Alice')
switch(partners, 'Bob', 'Alice')
switch(partners, 'Eve', 'Gertrude')
switch(partners, 'Fred', 'Claire')
switch(partners, 'Fred', 'Bob')
switch(partners, 'Dave', 'Fred')

# Output Fred's partner at the end
print(partners['Fred'])