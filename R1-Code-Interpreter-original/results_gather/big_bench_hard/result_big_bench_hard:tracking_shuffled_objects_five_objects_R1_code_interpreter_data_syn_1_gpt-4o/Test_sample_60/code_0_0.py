# Initial pairings
partners = {
    'Alice': 'Ophelia',
    'Bob': 'Karl',
    'Claire': 'Lola',
    'Dave': 'Jamie',
    'Eve': 'Izzi'
}

# Function to switch partners
def switch_partners(person1, person2):
    partners[person1], partners[person2] = partners[person2], partners[person1]

# Apply the partner switches
switch_partners('Claire', 'Alice')
switch_partners('Eve', 'Alice')
switch_partners('Claire', 'Dave')
switch_partners('Dave', 'Bob')
switch_partners('Claire', 'Alice')

# Output the final partner of Eve
print(partners['Eve'])