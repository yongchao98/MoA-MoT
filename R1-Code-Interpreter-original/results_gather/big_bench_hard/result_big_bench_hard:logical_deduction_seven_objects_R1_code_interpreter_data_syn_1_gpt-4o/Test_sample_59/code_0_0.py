from itertools import permutations

books = ['R', 'W', 'P', 'B', 'G', 'O', 'Bl']

# Conditions
def is_valid_order(order):
    return (order.index('R') == 5 and  # Red book is second from the right
            order.index('B') == 3 and  # Black book is fourth from the left
            order.index('R') > order.index('W') and  # Red is to the right of White
            order.index('G') < order.index('O') and  # Gray is to the left of Orange
            order.index('G') > order.index('Bl') and  # Gray is to the right of Blue
            order.index('O') < order.index('W'))  # Orange is to the left of White

# Find the valid order
for order in permutations(books):
    if is_valid_order(order):
        print(order)
        break