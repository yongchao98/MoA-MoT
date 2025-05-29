from itertools import permutations

# Define the books
books = ['orange', 'red', 'green', 'white', 'purple', 'brown', 'yellow']

# Define the constraints
def is_valid_order(order):
    green_index = order.index('green')
    purple_index = order.index('purple')
    white_index = order.index('white')
    orange_index = order.index('orange')
    brown_index = order.index('brown')
    yellow_index = order.index('yellow')
    red_index = order.index('red')
    
    return (green_index < purple_index and
            white_index > orange_index and
            green_index > brown_index and
            yellow_index > white_index and
            red_index > yellow_index and
            brown_index > red_index)

# Find the valid order
for order in permutations(books):
    if is_valid_order(order):
        print(order)
        break