# Initial pairings
partners = {
    'Alice': 'Lola',
    'Bob': 'Rodrigo',
    'Claire': 'Jamie',
    'Dave': 'Izzi',
    'Eve': 'Sam',
    'Fred': 'Karl',
    'Gertrude': 'Helga'
}

# Function to switch partners
def switch(partners, dancer1, dancer2):
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Apply the switches
switch(partners, 'Alice', 'Claire')
switch(partners, 'Bob', 'Fred')
switch(partners, 'Claire', 'Dave')
switch(partners, 'Alice', 'Fred')
switch(partners, 'Claire', 'Alice')
switch(partners, 'Bob', 'Eve')
switch(partners, 'Eve', 'Gertrude')

# Output the final partner of Eve
print(partners['Eve'])