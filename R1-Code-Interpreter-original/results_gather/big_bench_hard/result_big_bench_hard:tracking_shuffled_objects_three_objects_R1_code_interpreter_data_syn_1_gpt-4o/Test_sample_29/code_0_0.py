# Initial partners
partners = {
    'Alice': 'Karl',
    'Bob': 'Rodrigo',
    'Claire': 'Patrick'
}

# First switch: Claire and Alice switch partners
partners['Alice'], partners['Claire'] = partners['Claire'], partners['Alice']

# Second switch: Bob and Alice switch partners
partners['Alice'], partners['Bob'] = partners['Bob'], partners['Alice']

# Third switch: Alice and Claire switch partners
partners['Alice'], partners['Claire'] = partners['Claire'], partners['Alice']

# Output the final partner of Claire
print(partners['Claire'])