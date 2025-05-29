from itertools import product

# Define the numbers and target
numbers = [3, 1, 9, 6]
target = 24
operations = ['+', '-', '*', '/']

# Function to evaluate expression with given operations
def evaluate_expression(numbers, ops):
    try:
        # Construct the expression with operations
        expression = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}"
        # Evaluate the expression
        return eval(expression)
    except ZeroDivisionError:
        return None

# Try all combinations of operations
for ops in product(operations, repeat=3):
    # Evaluate without parentheses
    if evaluate_expression(numbers, ops) == target:
        print([numbers[0], ops[0], numbers[1], ops[1], numbers[2], ops[2], numbers[3]])
        break
    # Evaluate with different parentheses
    if evaluate_expression([evaluate_expression([numbers[0], numbers[1]], ops[:1]), numbers[2], numbers[3]], ops[1:]) == target:
        print(['(', numbers[0], ops[0], numbers[1], ')', ops[1], numbers[2], ops[2], numbers[3]])
        break
    if evaluate_expression([numbers[0], evaluate_expression([numbers[1], numbers[2]], ops[1:2]), numbers[3]], [ops[0], ops[2]]) == target:
        print([numbers[0], ops[0], '(', numbers[1], ops[1], numbers[2], ')', ops[2], numbers[3]])
        break
    if evaluate_expression([numbers[0], numbers[1], evaluate_expression([numbers[2], numbers[3]], ops[2:])], ops[:2]) == target:
        print([numbers[0], ops[0], numbers[1], ops[1], '(', numbers[2], ops[2], numbers[3], ')'])
        break