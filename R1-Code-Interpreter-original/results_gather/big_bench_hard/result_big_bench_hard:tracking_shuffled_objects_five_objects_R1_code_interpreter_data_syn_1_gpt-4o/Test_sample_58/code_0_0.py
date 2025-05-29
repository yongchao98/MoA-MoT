# Initial pairings
partners = {
    'Alice': 'Rodrigo',
    'Bob': 'Sam',
    'Claire': 'Karl',
    'Dave': 'Izzi',
    'Eve': 'Jamie'
}

# Switches
switches = [
    ('Eve', 'Bob'),
    ('Dave', 'Alice'),
    ('Alice', 'Claire'),
    ('Eve', 'Claire'),
    ('Alice', 'Bob')
]

# Apply switches
for dancer1, dancer2 in switches:
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Find Eve's final partner
final_partner_eve = partners['Eve']
print(final_partner_eve)