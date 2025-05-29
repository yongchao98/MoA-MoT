# Initial distribution of presents
presents = {
    'Alice': 'Green',
    'Bob': 'Brown',
    'Claire': 'Purple',
    'Dave': 'White',
    'Eve': 'Blue'
}

# Swap operations
swaps = [
    ('Alice', 'Bob'),
    ('Bob', 'Eve'),
    ('Dave', 'Eve'),
    ('Alice', 'Claire'),
    ('Bob', 'Eve')
]

# Perform the swaps
for person1, person2 in swaps:
    presents[person1], presents[person2] = presents[person2], presents[person1]

# Output the present Dave has at the end
print(presents['Dave'])