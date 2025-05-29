# Initial pairings
partners = {
    'Alice': 'Ophelia',
    'Bob': 'Izzi',
    'Claire': 'Helga',
    'Dave': 'Patrick',
    'Eve': 'Lola'
}

# Switches
switches = [
    ('Eve', 'Claire'),
    ('Claire', 'Dave'),
    ('Claire', 'Alice'),
    ('Bob', 'Alice'),
    ('Dave', 'Claire')
]

# Apply each switch
for dancer1, dancer2 in switches:
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Find Dave's partner
daves_partner = partners['Dave']
print(daves_partner)