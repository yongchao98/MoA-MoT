# Initial pairings
partners = {
    'Alice': 'Melissa',
    'Bob': 'Lola',
    'Claire': 'Patrick',
    'Dave': 'Sam',
    'Eve': 'Izzi',
    'Fred': 'Helga',
    'Gertrude': 'Rodrigo'
}

# Function to switch partners
def switch(partners, dancer1, dancer2):
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Apply the switches
switch(partners, 'Fred', 'Gertrude')
switch(partners, 'Claire', 'Fred')
switch(partners, 'Dave', 'Alice')
switch(partners, 'Alice', 'Bob')
switch(partners, 'Gertrude', 'Eve')
switch(partners, 'Dave', 'Gertrude')
switch(partners, 'Alice', 'Dave')

# Output the final partner of Eve
print(partners['Eve'])