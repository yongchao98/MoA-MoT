# Initial pairings
partners = {
    'Alice': 'Melissa',
    'Bob': 'Ophelia',
    'Claire': 'Patrick',
    'Dave': 'Jamie',
    'Eve': 'Lola',
    'Fred': 'Helga',
    'Gertrude': 'Izzi'
}

# List of switches
switches = [
    ('Alice', 'Fred'),
    ('Gertrude', 'Eve'),
    ('Gertrude', 'Alice'),
    ('Bob', 'Alice'),
    ('Dave', 'Fred'),
    ('Eve', 'Alice'),
    ('Eve', 'Claire')
]

# Apply each switch
for dancer1, dancer2 in switches:
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Output the final partner of Gertrude
print(partners['Gertrude'])