# Initial pairings
partners = {
    'Alice': 'Jamie',
    'Bob': 'Patrick',
    'Claire': 'Karl',
    'Dave': 'Helga',
    'Eve': 'Izzi'
}

# Switches
switches = [
    ('Claire', 'Alice'),
    ('Alice', 'Dave'),
    ('Alice', 'Bob'),
    ('Bob', 'Eve'),
    ('Dave', 'Bob')
]

# Perform the switches
for dancer1, dancer2 in switches:
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Find Eve's partner
eves_partner = partners['Eve']
print(eves_partner)