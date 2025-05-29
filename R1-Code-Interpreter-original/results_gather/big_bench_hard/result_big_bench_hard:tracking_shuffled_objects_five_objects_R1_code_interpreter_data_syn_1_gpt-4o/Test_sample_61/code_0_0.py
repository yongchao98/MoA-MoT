# Initial pairings
partners = {
    'Alice': 'Karl',
    'Bob': 'Patrick',
    'Claire': 'Jamie',
    'Dave': 'Rodrigo',
    'Eve': 'Lola'
}

# Partner switches
switches = [
    ('Eve', 'Claire'),
    ('Bob', 'Alice'),
    ('Dave', 'Bob'),
    ('Dave', 'Eve'),
    ('Claire', 'Eve')
]

# Apply each switch
for dancer1, dancer2 in switches:
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Find Dave's partner at the end
daves_partner = partners['Dave']
print(daves_partner)