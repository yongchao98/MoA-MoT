# Initial partners
partners = {
    'Alice': 'Lola',
    'Bob': 'Patrick',
    'Claire': 'Melissa'
}

# First switch: Alice and Claire switch partners
partners['Alice'], partners['Claire'] = partners['Claire'], partners['Alice']

# Second switch: Bob and Claire switch partners
partners['Bob'], partners['Claire'] = partners['Claire'], partners['Bob']

# Third switch: Claire and Alice switch partners
partners['Claire'], partners['Alice'] = partners['Alice'], partners['Claire']

# Output Bob's partner
print(partners['Bob'])