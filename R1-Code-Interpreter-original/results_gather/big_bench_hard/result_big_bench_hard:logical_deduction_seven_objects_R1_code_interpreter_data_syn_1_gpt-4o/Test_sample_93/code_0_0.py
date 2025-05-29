from itertools import permutations

books = ['green', 'brown', 'white', 'black', 'orange', 'purple', 'yellow']

# Conditions
def is_valid_order(order):
    return (order.index('white') == 3 and
            order.index('black') == 1 and
            order.index('purple') < order.index('yellow') and
            order.index('green') < order.index('white') and
            order.index('brown') > order.index('yellow') and
            order.index('green') > order.index('orange'))

# Find the valid order
for order in permutations(books):
    if is_valid_order(order):
        print(order)
        break