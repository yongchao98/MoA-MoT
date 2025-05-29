# Initial presents
presents = {
    'Alice': 'Green',
    'Bob': 'Yellow',
    'Claire': 'Orange',
    'Dave': 'White',
    'Eve': 'Brown'
}

# Swap operations
swaps = [
    ('Eve', 'Dave'),
    ('Dave', 'Alice'),
    ('Bob', 'Dave'),
    ('Eve', 'Claire'),
    ('Eve', 'Alice')
]

# Perform swaps
for person1, person2 in swaps:
    presents[person1], presents[person2] = presents[person2], presents[person1]

# Output the present Dave has at the end
print(presents['Dave'])