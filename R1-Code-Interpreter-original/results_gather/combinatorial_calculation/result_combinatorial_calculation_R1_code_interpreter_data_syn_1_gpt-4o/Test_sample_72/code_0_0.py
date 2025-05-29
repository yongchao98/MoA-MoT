from itertools import product

# Given numbers
numbers = [4, 4, 7, 8]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all combinations of operations and parentheses
for ops in product(operations, repeat=3):
    # Without parentheses
    expr1 = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}"
    if evaluate_expression(expr1) == 24:
        print([numbers[0], ops[0], numbers[1], ops[1], numbers[2], ops[2], numbers[3]])
        break
    
    # With parentheses around first two numbers
    expr2 = f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}{numbers[2]}{ops[2]}{numbers[3]}"
    if evaluate_expression(expr2) == 24:
        print(['(', numbers[0], ops[0], numbers[1], ')', ops[1], numbers[2], ops[2], numbers[3]])
        break
    
    # With parentheses around last two numbers
    expr3 = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}({numbers[2]}{ops[2]}{numbers[3]})"
    if evaluate_expression(expr3) == 24:
        print([numbers[0], ops[0], numbers[1], ops[1], '(', numbers[2], ops[2], numbers[3], ')'])
        break
    
    # With parentheses around middle two numbers
    expr4 = f"{numbers[0]}{ops[0]}({numbers[1]}{ops[1]}{numbers[2]}){ops[2]}{numbers[3]}"
    if evaluate_expression(expr4) == 24:
        print([numbers[0], ops[0], '(', numbers[1], ops[1], numbers[2], ')', ops[2], numbers[3]])
        break
    
    # With parentheses around first three numbers
    expr5 = f"({numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}){ops[2]}{numbers[3]}"
    if evaluate_expression(expr5) == 24:
        print(['(', numbers[0], ops[0], numbers[1], ops[1], numbers[2], ')', ops[2], numbers[3]])
        break
    
    # With parentheses around last three numbers
    expr6 = f"{numbers[0]}{ops[0]}({numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]})"
    if evaluate_expression(expr6) == 24:
        print([numbers[0], ops[0], '(', numbers[1], ops[1], numbers[2], ops[2], numbers[3], ')'])
        break