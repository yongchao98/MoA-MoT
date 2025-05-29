# Initial partners
partners = {
    'Alice': 'Lola',
    'Bob': 'Melissa',
    'Claire': 'Sam',
    'Dave': 'Karl',
    'Eve': 'Helga'
}

# Switches
switches = [
    ('Dave', 'Alice'),
    ('Claire', 'Dave'),
    ('Alice', 'Dave'),
    ('Bob', 'Eve'),
    ('Bob', 'Alice')
]

# Apply switches
for dancer1, dancer2 in switches:
    partners[dancer1], partners[dancer2] = partners[dancer2], partners[dancer1]

# Final partner of Bob
print(partners['Bob'])