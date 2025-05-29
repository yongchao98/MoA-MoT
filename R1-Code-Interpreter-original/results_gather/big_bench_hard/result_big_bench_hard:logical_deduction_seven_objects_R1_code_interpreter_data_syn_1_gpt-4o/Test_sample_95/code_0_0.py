from itertools import permutations

fruits = ['apples', 'plums', 'mangoes', 'watermelons', 'oranges', 'kiwis', 'loquats']

# Constraints
def is_valid_order(order):
    mangoes_index = order.index('mangoes')
    plums_index = order.index('plums')
    watermelons_index = order.index('watermelons')
    oranges_index = order.index('oranges')
    kiwis_index = order.index('kiwis')
    apples_index = order.index('apples')
    loquats_index = order.index('loquats')
    
    return (mangoes_index == 2 and
            plums_index > mangoes_index and
            watermelons_index < plums_index and
            watermelons_index > oranges_index and
            kiwis_index < oranges_index and
            apples_index > loquats_index)

# Find the valid order
for order in permutations(fruits):
    if is_valid_order(order):
        print(order)
        break