# Initial pairings
partners = {
    'Alice': 'Karl',
    'Bob': 'Helga',
    'Claire': 'Sam',
    'Dave': 'Rodrigo',
    'Eve': 'Melissa',
    'Fred': 'Izzi',
    'Gertrude': 'Ophelia'
}

# Function to switch partners
def switch(partners, person1, person2):
    partners[person1], partners[person2] = partners[person2], partners[person1]

# Apply the switches
switch(partners, 'Gertrude', 'Fred')
switch(partners, 'Gertrude', 'Bob')
switch(partners, 'Alice', 'Bob')
switch(partners, 'Eve', 'Dave')
switch(partners, 'Bob', 'Gertrude')
switch(partners, 'Alice', 'Claire')
switch(partners, 'Dave', 'Fred')

# Find Bob's final partner
final_partner_of_bob = partners['Bob']
print(final_partner_of_bob)