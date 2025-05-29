# Initial pairings
partners = {
    'Alice': 'Jamie',
    'Bob': 'Lola',
    'Claire': 'Izzi',
    'Dave': 'Rodrigo',
    'Eve': 'Ophelia'
}

# Switches
switches = [
    ('Bob', 'Eve'),
    ('Alice', 'Bob'),
    ('Dave', 'Alice'),
    ('Dave', 'Claire'),
    ('Bob', 'Claire')
]

# Apply switches
for dancer1, dancer2 in switches:
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Find Claire's partner
claire_partner = partners['Claire']
print(claire_partner)