# Initial pairings
partners = {
    'Alice': 'Karl',
    'Bob': 'Helga',
    'Claire': 'Sam',
    'Dave': 'Rodrigo',
    'Eve': 'Melissa',
    'Fred': 'Izzi',
    'Gertrude': 'Ophelia'
}

# Function to switch partners
def switch(partners, dancer1, dancer2):
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Apply the switches
switch(partners, 'Gertrude', 'Fred')
switch(partners, 'Gertrude', 'Bob')
switch(partners, 'Alice', 'Bob')
switch(partners, 'Eve', 'Dave')
switch(partners, 'Bob', 'Gertrude')
switch(partners, 'Alice', 'Claire')
switch(partners, 'Dave', 'Fred')

# Output who Bob is dancing with
print(partners['Bob'])