from itertools import product

# Define the numbers and target
numbers = [6, 2, 7]
target = 10

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all possible expressions
for ops in product(operations, repeat=2):
    # Without parentheses
    expr1 = f"{numbers[0]} {ops[0]} {numbers[1]} {ops[1]} {numbers[2]}"
    if evaluate_expression(expr1) == target:
        print([numbers[0], ops[0], numbers[1], ops[1], numbers[2]])
    
    # With parentheses
    expr2 = f"({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]}"
    if evaluate_expression(expr2) == target:
        print(['(', numbers[0], ops[0], numbers[1], ')', ops[1], numbers[2]])
    
    expr3 = f"{numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} {numbers[2]})"
    if evaluate_expression(expr3) == target:
        print([numbers[0], ops[0], '(', numbers[1], ops[1], numbers[2], ')'])