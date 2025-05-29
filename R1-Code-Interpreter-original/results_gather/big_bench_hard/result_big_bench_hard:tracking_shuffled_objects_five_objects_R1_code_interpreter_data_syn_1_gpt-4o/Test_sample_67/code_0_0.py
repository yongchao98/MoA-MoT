# Initial pairings
partners = {
    'Alice': 'Karl',
    'Bob': 'Patrick',
    'Claire': 'Jamie',
    'Dave': 'Rodrigo',
    'Eve': 'Lola'
}

# Function to switch partners
def switch_partners(person1, person2):
    partners[person1], partners[person2] = partners[person2], partners[person1]

# Apply the partner switches
switch_partners('Eve', 'Claire')
switch_partners('Bob', 'Alice')
switch_partners('Dave', 'Bob')
switch_partners('Dave', 'Eve')
switch_partners('Claire', 'Eve')

# Find out who Eve is dancing with
final_partner_of_eve = partners['Eve']
print(final_partner_of_eve)