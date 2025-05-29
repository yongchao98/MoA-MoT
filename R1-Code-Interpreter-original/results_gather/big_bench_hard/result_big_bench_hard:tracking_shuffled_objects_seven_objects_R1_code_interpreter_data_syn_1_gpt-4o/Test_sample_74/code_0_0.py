# Initial pairings
partners = {
    'Alice': 'Karl',
    'Bob': 'Izzi',
    'Claire': 'Melissa',
    'Dave': 'Lola',
    'Eve': 'Rodrigo',
    'Fred': 'Helga',
    'Gertrude': 'Patrick'
}

# List of partner switches
switches = [
    ('Alice', 'Bob'),
    ('Eve', 'Gertrude'),
    ('Claire', 'Gertrude'),
    ('Gertrude', 'Bob'),
    ('Claire', 'Bob'),
    ('Fred', 'Gertrude'),
    ('Claire', 'Dave')
]

# Apply each switch
for dancer1, dancer2 in switches:
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Output the final partner of Bob
print(partners['Bob'])