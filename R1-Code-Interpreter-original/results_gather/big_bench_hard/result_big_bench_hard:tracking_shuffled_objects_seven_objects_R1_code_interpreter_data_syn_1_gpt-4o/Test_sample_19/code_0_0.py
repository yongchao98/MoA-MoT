# Initial pairings
partners = {
    'Alice': 'Karl',
    'Bob': 'Izzi',
    'Claire': 'Jamie',
    'Dave': 'Melissa',
    'Eve': 'Lola',
    'Fred': 'Patrick',
    'Gertrude': 'Rodrigo'
}

# Function to switch partners
def switch(partners, dancer1, dancer2):
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Apply the switches
switch(partners, 'Alice', 'Eve')
switch(partners, 'Eve', 'Bob')
switch(partners, 'Gertrude', 'Dave')
switch(partners, 'Eve', 'Claire')
switch(partners, 'Alice', 'Bob')
switch(partners, 'Bob', 'Fred')
switch(partners, 'Fred', 'Claire')

# Output Fred's partner at the end
print(partners['Fred'])