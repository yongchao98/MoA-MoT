# Initial pairings
partners = {
    'Alice': 'Patrick',
    'Bob': 'Rodrigo',
    'Claire': 'Lola',
    'Dave': 'Izzi',
    'Eve': 'Sam',
    'Fred': 'Karl',
    'Gertrude': 'Melissa'
}

# Function to switch partners
def switch(partners, dancer1, dancer2):
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Apply the switches
switch(partners, 'Bob', 'Fred')
switch(partners, 'Claire', 'Alice')
switch(partners, 'Claire', 'Dave')
switch(partners, 'Gertrude', 'Alice')
switch(partners, 'Claire', 'Bob')
switch(partners, 'Dave', 'Fred')
switch(partners, 'Eve', 'Alice')

# Output the partner of Eve
print(partners['Eve'])