# Initial distribution of presents
presents = {
    'Alice': 'Brown',
    'Bob': 'White',
    'Claire': 'Black',
    'Dave': 'Red',
    'Eve': 'Purple'
}

# Swap 1: Bob and Claire
presents['Bob'], presents['Claire'] = presents['Claire'], presents['Bob']

# Swap 2: Eve and Alice
presents['Eve'], presents['Alice'] = presents['Alice'], presents['Eve']

# Swap 3: Bob and Claire
presents['Bob'], presents['Claire'] = presents['Claire'], presents['Bob']

# Swap 4: Claire and Dave
presents['Claire'], presents['Dave'] = presents['Dave'], presents['Claire']

# Swap 5: Eve and Dave
presents['Eve'], presents['Dave'] = presents['Dave'], presents['Eve']

# Output the present Eve has at the end
print(presents['Eve'])