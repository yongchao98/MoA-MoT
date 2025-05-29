from itertools import product

# Given numbers
numbers = [1, 9, 4, 9, 7]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression with parentheses
def evaluate_expression(numbers, ops, paren_pos):
    expr = []
    for i, num in enumerate(numbers):
        if i in paren_pos:
            expr.append('(')
        expr.append(str(num))
        if i in paren_pos:
            expr.append(')')
        if i < len(ops):
            expr.append(ops[i])
    expr_str = ''.join(expr)
    try:
        return eval(expr_str)
    except ZeroDivisionError:
        return None

# Try all combinations of operations and parentheses
for ops in product(operations, repeat=4):
    for paren_pos in product(range(5), repeat=2):
        if paren_pos[0] < paren_pos[1]:  # Ensure valid parentheses
            result = evaluate_expression(numbers, ops, paren_pos)
            if result == 69:
                expr = []
                for i, num in enumerate(numbers):
                    if i in paren_pos:
                        expr.append('(')
                    expr.append(num)
                    if i in paren_pos:
                        expr.append(')')
                    if i < len(ops):
                        expr.append(ops[i])
                print(expr)
                break