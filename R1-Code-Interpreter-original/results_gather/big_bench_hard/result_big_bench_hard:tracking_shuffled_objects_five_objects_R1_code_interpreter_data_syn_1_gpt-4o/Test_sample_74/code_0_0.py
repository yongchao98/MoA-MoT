# Initial partners
partners = {
    'Alice': 'Sam',
    'Bob': 'Izzi',
    'Claire': 'Karl',
    'Dave': 'Helga',
    'Eve': 'Ophelia'
}

# Function to switch partners
def switch_partners(person1, person2):
    partners[person1], partners[person2] = partners[person2], partners[person1]

# Apply the partner switches
switch_partners('Bob', 'Alice')
switch_partners('Alice', 'Claire')
switch_partners('Alice', 'Dave')
switch_partners('Alice', 'Eve')
switch_partners('Dave', 'Claire')

# Output the final partner of Claire
print(partners['Claire'])