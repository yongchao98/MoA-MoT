# Initial partners
partners = {
    'Alice': 'Helga',
    'Bob': 'Ophelia',
    'Claire': 'Sam'
}

# First switch: Bob and Alice switch partners
partners['Alice'], partners['Bob'] = partners['Bob'], partners['Alice']

# Second switch: Claire and Bob switch partners
partners['Claire'], partners['Bob'] = partners['Bob'], partners['Claire']

# Third switch: Claire and Alice switch partners
partners['Claire'], partners['Alice'] = partners['Alice'], partners['Claire']

# Final partner of Claire
print(partners['Claire'])