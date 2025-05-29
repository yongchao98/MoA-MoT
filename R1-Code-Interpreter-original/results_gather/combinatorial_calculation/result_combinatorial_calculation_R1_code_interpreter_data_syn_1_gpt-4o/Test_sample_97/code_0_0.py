from itertools import product

# The numbers in the given order
numbers = [9, 9, 1, 1, 1]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression given as a list
def evaluate_expression(expr):
    try:
        return eval(''.join(map(str, expr)))
    except ZeroDivisionError:
        return None

# Try all combinations of operations and parentheses
for ops in product(operations, repeat=4):
    # Generate possible expressions with parentheses
    expressions = [
        [numbers[0], ops[0], numbers[1], ops[1], numbers[2], ops[2], numbers[3], ops[3], numbers[4]],
        ['(', numbers[0], ops[0], numbers[1], ')', ops[1], numbers[2], ops[2], numbers[3], ops[3], numbers[4]],
        [numbers[0], ops[0], '(', numbers[1], ops[1], numbers[2], ')', ops[2], numbers[3], ops[3], numbers[4]],
        [numbers[0], ops[0], numbers[1], ops[1], '(', numbers[2], ops[2], numbers[3], ')', ops[3], numbers[4]],
        [numbers[0], ops[0], numbers[1], ops[1], numbers[2], ops[2], '(', numbers[3], ops[3], numbers[4], ')'],
        ['(', numbers[0], ops[0], numbers[1], ops[1], numbers[2], ')', ops[2], numbers[3], ops[3], numbers[4]],
        [numbers[0], ops[0], '(', numbers[1], ops[1], numbers[2], ops[2], numbers[3], ')', ops[3], numbers[4]],
        ['(', numbers[0], ops[0], numbers[1], ')', ops[1], '(', numbers[2], ops[2], numbers[3], ')', ops[3], numbers[4]],
    ]
    
    for expr in expressions:
        if evaluate_expression(expr) == 78:
            print(expr)
            break