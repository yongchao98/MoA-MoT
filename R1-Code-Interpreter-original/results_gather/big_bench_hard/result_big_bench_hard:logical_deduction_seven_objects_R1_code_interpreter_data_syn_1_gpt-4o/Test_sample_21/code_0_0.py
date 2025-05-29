from itertools import permutations

birds = ['hawk', 'hummingbird', 'quail', 'owl', 'crow', 'robin', 'cardinal']

# Constraints
def is_valid_order(order):
    return (order[6] == 'hawk' and  # hawk is the rightmost
            order[5] == 'cardinal' and  # cardinal is second from the right
            order[3] == 'robin' and  # robin is fourth from the left
            order.index('owl') < order.index('hummingbird') and  # owl is to the left of hummingbird
            order.index('robin') < order.index('quail') and  # robin is to the left of quail
            order.index('owl') > order.index('crow'))  # owl is to the right of crow

# Find the valid order
for order in permutations(birds):
    if is_valid_order(order):
        print(order)
        break