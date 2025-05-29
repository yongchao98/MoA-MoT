# Initial distribution of presents
presents = {
    'Alice': 'Green',
    'Bob': 'Yellow',
    'Claire': 'Orange',
    'Dave': 'White',
    'Eve': 'Brown'
}

# Swap 1: Eve and Dave
presents['Eve'], presents['Dave'] = presents['Dave'], presents['Eve']

# Swap 2: Dave and Alice
presents['Dave'], presents['Alice'] = presents['Alice'], presents['Dave']

# Swap 3: Bob and Dave
presents['Bob'], presents['Dave'] = presents['Dave'], presents['Bob']

# Swap 4: Eve and Claire
presents['Eve'], presents['Claire'] = presents['Claire'], presents['Eve']

# Swap 5: Eve and Alice
presents['Eve'], presents['Alice'] = presents['Alice'], presents['Eve']

# Output the present Bob has at the end
print(presents['Bob'])