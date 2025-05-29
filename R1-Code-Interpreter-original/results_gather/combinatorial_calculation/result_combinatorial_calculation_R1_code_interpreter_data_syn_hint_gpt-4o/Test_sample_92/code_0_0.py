from itertools import product

# Given numbers
numbers = [1, 8, 5, 2, 3]
# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression with given operations and parentheses
def evaluate_expression(ops):
    try:
        # Construct the expression with operations
        expression = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}"
        # Evaluate the expression
        return eval(expression)
    except ZeroDivisionError:
        return None

# Try all combinations of operations
for ops in product(operations, repeat=4):
    # Check if the expression evaluates to 78
    if evaluate_expression(ops) == 78:
        # Print the successful combination
        print([numbers[0], ops[0], numbers[1], ops[1], numbers[2], ops[2], numbers[3], ops[3], numbers[4]])
        break