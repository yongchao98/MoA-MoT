# Initial partners
partners = {
    'Alice': 'Helga',
    'Bob': 'Patrick',
    'Claire': 'Sam'
}

# First switch: Claire and Alice switch partners
partners['Alice'], partners['Claire'] = partners['Claire'], partners['Alice']

# Second switch: Claire and Bob switch partners
partners['Claire'], partners['Bob'] = partners['Bob'], partners['Claire']

# Third switch: Alice and Bob switch partners
partners['Alice'], partners['Bob'] = partners['Bob'], partners['Alice']

# Final partner of Alice
print(partners['Alice'])