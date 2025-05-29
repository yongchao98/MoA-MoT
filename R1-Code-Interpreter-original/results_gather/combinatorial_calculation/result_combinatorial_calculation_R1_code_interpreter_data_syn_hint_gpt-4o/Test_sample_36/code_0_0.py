from itertools import product

# Given numbers
numbers = [9, 6, 5, 7, 7]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression with parentheses
def evaluate_expression(numbers, ops, parens):
    expr = f"{numbers[0]}"
    for i in range(4):
        expr += f" {ops[i]} {numbers[i+1]}"
    if parens:
        expr = f"({expr})"
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all combinations of operations and parentheses
for ops in product(operations, repeat=4):
    for parens in [True, False]:
        result = evaluate_expression(numbers, ops, parens)
        if result == 69:
            expr = f"{numbers[0]}"
            for i in range(4):
                expr += f" {ops[i]} {numbers[i+1]}"
            if parens:
                expr = f"({expr})"
            print(expr)
            break