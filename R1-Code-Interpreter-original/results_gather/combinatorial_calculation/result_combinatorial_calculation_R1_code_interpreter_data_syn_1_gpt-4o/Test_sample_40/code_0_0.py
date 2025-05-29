from itertools import product

# Given numbers
numbers = [9, 5, 6, 4, 6, 9]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression with given operations and parentheses
def evaluate_expression(ops):
    try:
        # Construct the expression with operations
        expression = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}{ops[4]}{numbers[5]}"
        # Evaluate the expression
        return eval(expression)
    except ZeroDivisionError:
        return None

# Try all combinations of operations
for ops in product(operations, repeat=5):
    if evaluate_expression(ops) == 100:
        # If the expression evaluates to 100, print the operations
        print([numbers[0], ops[0], numbers[1], ops[1], numbers[2], ops[2], numbers[3], ops[3], numbers[4], ops[4], numbers[5]])
        break