from itertools import product

# Given numbers
numbers = [8, 5, 3, 2, 5]

# Possible operations
operations = ['+', '-', '*', '/']

# Target value
target = 69

# Function to evaluate expression with given operations and parentheses
def evaluate_expression(numbers, ops, parens):
    expr = f"{parens[0]}{numbers[0]}{ops[0]}{parens[1]}{numbers[1]}{ops[1]}{parens[2]}{numbers[2]}{ops[2]}{parens[3]}{numbers[3]}{ops[3]}{numbers[4]}{parens[4]}"
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all combinations of operations and parentheses
for ops in product(operations, repeat=4):
    for parens in product(['', '(', ')'], repeat=5):
        if parens.count('(') == parens.count(')'):
            result = evaluate_expression(numbers, ops, parens)
            if result == target:
                print([numbers[0], ops[0], parens[1], numbers[1], ops[1], parens[2], numbers[2], ops[2], parens[3], numbers[3], ops[3], numbers[4], parens[4]])
                break