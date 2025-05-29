# Initial pairings
partners = {
    'Alice': 'Lola',
    'Bob': 'Melissa',
    'Claire': 'Helga',
    'Dave': 'Karl',
    'Eve': 'Sam',
    'Fred': 'Izzi',
    'Gertrude': 'Patrick'
}

# Function to switch partners
def switch(partners, dancer1, dancer2):
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Apply the switches
switch(partners, 'Alice', 'Eve')
switch(partners, 'Dave', 'Fred')
switch(partners, 'Eve', 'Claire')
switch(partners, 'Dave', 'Gertrude')
switch(partners, 'Dave', 'Bob')
switch(partners, 'Alice', 'Claire')
switch(partners, 'Eve', 'Gertrude')

# Output the partner of Bob
print(partners['Bob'])