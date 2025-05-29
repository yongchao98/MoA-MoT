# Initial partners
partners = {
    'Alice': 'Sam',
    'Bob': 'Lola',
    'Claire': 'Karl'
}

# First switch: Bob and Alice switch partners
partners['Alice'], partners['Bob'] = partners['Bob'], partners['Alice']

# Second switch: Bob and Claire switch partners
partners['Bob'], partners['Claire'] = partners['Claire'], partners['Bob']

# Third switch: Alice and Bob switch partners
partners['Alice'], partners['Bob'] = partners['Bob'], partners['Alice']

# Final partner of Claire
print(partners['Claire'])