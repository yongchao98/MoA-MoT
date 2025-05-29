from itertools import product

# Given numbers
numbers = [5, 5, 9]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all combinations of operations and parentheses
for ops in product(operations, repeat=2):
    # Without parentheses
    expr1 = f"{numbers[0]} {ops[0]} {numbers[1]} {ops[1]} {numbers[2]}"
    if evaluate_expression(expr1) == 10:
        print([numbers[0], ops[0], numbers[1], ops[1], numbers[2]])
    
    # With parentheses around first two numbers
    expr2 = f"({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]}"
    if evaluate_expression(expr2) == 10:
        print(['(', numbers[0], ops[0], numbers[1], ')', ops[1], numbers[2]])
    
    # With parentheses around last two numbers
    expr3 = f"{numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} {numbers[2]})"
    if evaluate_expression(expr3) == 10:
        print([numbers[0], ops[0], '(', numbers[1], ops[1], numbers[2], ')'])