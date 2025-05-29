# Initial pairings
partners = {
    'Alice': 'Patrick',
    'Bob': 'Lola',
    'Claire': 'Izzi',
    'Dave': 'Rodrigo',
    'Eve': 'Helga',
    'Fred': 'Sam',
    'Gertrude': 'Melissa'
}

# Function to switch partners
def switch(partners, dancer1, dancer2):
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Apply the switches
switch(partners, 'Eve', 'Bob')
switch(partners, 'Eve', 'Dave')
switch(partners, 'Fred', 'Gertrude')
switch(partners, 'Gertrude', 'Alice')
switch(partners, 'Alice', 'Dave')
switch(partners, 'Claire', 'Alice')
switch(partners, 'Alice', 'Dave')

# Output the final partner of Eve
print(partners['Eve'])