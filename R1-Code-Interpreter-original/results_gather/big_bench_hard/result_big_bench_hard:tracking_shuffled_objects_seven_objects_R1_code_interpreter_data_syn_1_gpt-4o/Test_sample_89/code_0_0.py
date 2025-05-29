# Initial pairings
partners = {
    'Alice': 'Patrick',
    'Bob': 'Ophelia',
    'Claire': 'Lola',
    'Dave': 'Rodrigo',
    'Eve': 'Jamie',
    'Fred': 'Sam',
    'Gertrude': 'Karl'
}

# List of partner switches
switches = [
    ('Claire', 'Fred'),
    ('Claire', 'Dave'),
    ('Alice', 'Bob'),
    ('Claire', 'Eve'),
    ('Dave', 'Claire'),
    ('Gertrude', 'Eve'),
    ('Bob', 'Eve')
]

# Function to switch partners
def switch_partners(partners, person1, person2):
    partners[person1], partners[person2] = partners[person2], partners[person1]

# Apply each switch
for person1, person2 in switches:
    switch_partners(partners, person1, person2)

# Output the final partner of Eve
print(partners['Eve'])