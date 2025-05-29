# Initial pairings
partners = {
    'Alice': 'Sam',
    'Bob': 'Rodrigo',
    'Claire': 'Jamie',
    'Dave': 'Lola',
    'Eve': 'Helga',
    'Fred': 'Melissa',
    'Gertrude': 'Ophelia'
}

# Function to switch partners
def switch(partners, dancer1, dancer2):
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Apply the partner switches
switch(partners, 'Bob', 'Dave')
switch(partners, 'Dave', 'Claire')
switch(partners, 'Claire', 'Eve')
switch(partners, 'Alice', 'Gertrude')
switch(partners, 'Alice', 'Bob')
switch(partners, 'Bob', 'Fred')
switch(partners, 'Fred', 'Eve')

# Output the final partner of Dave
print(partners['Dave'])