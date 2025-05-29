# Initial pairings
partners = {
    'Alice': 'Ophelia',
    'Bob': 'Rodrigo',
    'Claire': 'Izzi',
    'Dave': 'Melissa',
    'Eve': 'Lola',
    'Fred': 'Patrick',
    'Gertrude': 'Jamie'
}

# List of partner switches
switches = [
    ('Dave', 'Eve'),
    ('Gertrude', 'Fred'),
    ('Dave', 'Claire'),
    ('Gertrude', 'Alice'),
    ('Dave', 'Eve'),
    ('Bob', 'Alice'),
    ('Bob', 'Gertrude')
]

# Function to switch partners
def switch_partners(partners, person1, person2):
    partners[person1], partners[person2] = partners[person2], partners[person1]

# Apply each switch
for person1, person2 in switches:
    switch_partners(partners, person1, person2)

# Output the final partner of Eve
print(partners['Eve'])