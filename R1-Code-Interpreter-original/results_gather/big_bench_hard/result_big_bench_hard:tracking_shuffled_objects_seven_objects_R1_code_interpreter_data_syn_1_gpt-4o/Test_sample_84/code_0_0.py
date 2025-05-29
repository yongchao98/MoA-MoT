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
def switch_partners(dancer1, dancer2):
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Apply the partner switches
switch_partners('Eve', 'Bob')
switch_partners('Eve', 'Dave')
switch_partners('Fred', 'Gertrude')
switch_partners('Gertrude', 'Alice')
switch_partners('Alice', 'Dave')
switch_partners('Claire', 'Alice')
switch_partners('Alice', 'Dave')

# Output the final partner of Gertrude
print(partners['Gertrude'])