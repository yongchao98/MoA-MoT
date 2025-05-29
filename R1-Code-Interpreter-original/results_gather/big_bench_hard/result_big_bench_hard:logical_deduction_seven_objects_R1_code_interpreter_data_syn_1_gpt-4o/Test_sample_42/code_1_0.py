from itertools import permutations

# Define the fruits
fruits = ['mangoes', 'kiwis', 'oranges', 'watermelons', 'apples', 'pears', 'loquats']

# Define the constraints
def is_valid_order(order):
    # Apples are the cheapest
    if order.index('apples') != 6:
        return False
    # Oranges are the third-most expensive
    if order.index('oranges') != 2:
        return False
    # Loquats are less expensive than oranges
    if order.index('loquats') >= order.index('oranges'):
        return False
    # Kiwis are more expensive than pears
    if order.index('kiwis') <= order.index('pears'):
        return False
    # Mangoes are more expensive than watermelons
    if order.index('mangoes') <= order.index('watermelons'):
        return False
    # Loquats are more expensive than kiwis
    if order.index('loquats') <= order.index('kiwis'):
        return False
    return True

# Find the valid order
for order in permutations(fruits):
    if is_valid_order(order):
        print(order)
        break