from itertools import permutations

# Define the golfers
golfers = ['Mel', 'Rob', 'Eli', 'Dan', 'Ana', 'Ada', 'Mya']

# Define the constraints
def is_valid_order(order):
    return (order.index('Ana') == 6 and  # Ana is last
            order.index('Rob') == 3 and  # Rob is fourth
            order.index('Eli') == 1 and  # Eli is second
            order.index('Mya') < order.index('Rob') and  # Mya is above Rob
            order.index('Dan') < order.index('Mya') and  # Dan is above Mya
            order.index('Mel') < order.index('Ada'))  # Mel is above Ada

# Find the valid order
for order in permutations(golfers):
    if is_valid_order(order):
        print(order)
        break