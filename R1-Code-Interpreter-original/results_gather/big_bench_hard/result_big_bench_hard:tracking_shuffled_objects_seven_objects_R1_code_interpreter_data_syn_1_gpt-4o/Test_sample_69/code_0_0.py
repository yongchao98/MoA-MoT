# Initial pairings
partners = {
    'Alice': 'Lola',
    'Bob': 'Ophelia',
    'Claire': 'Sam',
    'Dave': 'Helga',
    'Eve': 'Patrick',
    'Fred': 'Jamie',
    'Gertrude': 'Rodrigo'
}

# Function to switch partners
def switch(partners, dancer1, dancer2):
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Perform the switches
switch(partners, 'Fred', 'Bob')
switch(partners, 'Eve', 'Alice')
switch(partners, 'Fred', 'Eve')
switch(partners, 'Claire', 'Fred')
switch(partners, 'Gertrude', 'Alice')
switch(partners, 'Alice', 'Bob')
switch(partners, 'Dave', 'Fred')

# Output the partner of Dave
print(partners['Dave'])