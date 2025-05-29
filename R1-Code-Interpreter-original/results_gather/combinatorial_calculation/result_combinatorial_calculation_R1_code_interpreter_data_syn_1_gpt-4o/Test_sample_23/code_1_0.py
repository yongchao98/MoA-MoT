from itertools import product

# Given numbers
numbers = [2, 2, 7, 2, 3]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression with given operations and parentheses
def evaluate_expression(ops):
    try:
        # Construct the expression with the current operations
        expression = f"(({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]}) {ops[2]} ({numbers[3]} {ops[3]} {numbers[4]})"
        # Evaluate the expression
        return eval(expression)
    except ZeroDivisionError:
        return None

# Iterate over all possible combinations of operations
for ops in product(operations, repeat=4):
    if evaluate_expression(ops) == 78:
        # If the expression evaluates to 78, print the operations
        print([numbers[0], ops[0], numbers[1], ops[1], numbers[2], ops[2], numbers[3], ops[3], numbers[4]])
        break