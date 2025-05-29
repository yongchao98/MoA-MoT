# Initial partners
partners = {
    'Alice': 'Izzi',
    'Bob': 'Melissa',
    'Claire': 'Lola'
}

# First switch: Claire and Alice switch partners
partners['Alice'], partners['Claire'] = partners['Claire'], partners['Alice']

# Second switch: Bob and Claire switch partners
partners['Bob'], partners['Claire'] = partners['Claire'], partners['Bob']

# Third switch: Alice and Bob switch partners
partners['Alice'], partners['Bob'] = partners['Bob'], partners['Alice']

# Final partner of Alice
print(partners['Alice'])