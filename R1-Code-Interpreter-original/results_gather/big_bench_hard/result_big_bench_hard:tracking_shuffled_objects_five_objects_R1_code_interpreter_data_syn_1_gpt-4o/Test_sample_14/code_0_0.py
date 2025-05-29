# Initial pairings
partners = {
    'Alice': 'Ophelia',
    'Bob': 'Izzi',
    'Claire': 'Lola',
    'Dave': 'Sam',
    'Eve': 'Melissa'
}

# Apply the partner switches
# 1. Claire and Bob switch partners
partners['Claire'], partners['Bob'] = partners['Bob'], partners['Claire']

# 2. Claire and Eve switch partners
partners['Claire'], partners['Eve'] = partners['Eve'], partners['Claire']

# 3. Bob and Dave switch partners
partners['Bob'], partners['Dave'] = partners['Dave'], partners['Bob']

# 4. Claire and Eve switch partners
partners['Claire'], partners['Eve'] = partners['Eve'], partners['Claire']

# 5. Alice and Eve switch partners
partners['Alice'], partners['Eve'] = partners['Eve'], partners['Alice']

# Determine who Dave is dancing with
print(partners['Dave'])