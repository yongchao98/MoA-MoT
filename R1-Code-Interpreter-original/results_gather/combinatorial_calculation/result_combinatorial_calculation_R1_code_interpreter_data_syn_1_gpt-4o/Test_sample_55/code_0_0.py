from itertools import product

# Given numbers
numbers = [3, 3, 9]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression with given operations
def evaluate_expression(nums, ops):
    try:
        # Construct the expression
        expression = f"{nums[0]}{ops[0]}{nums[1]}{ops[1]}{nums[2]}"
        # Evaluate the expression
        return eval(expression)
    except ZeroDivisionError:
        return None

# Try all combinations of operations
for ops in product(operations, repeat=2):
    # Evaluate without parentheses
    if evaluate_expression(numbers, ops) == 10:
        print([numbers[0], ops[0], numbers[1], ops[1], numbers[2]])
    
    # Evaluate with parentheses around first two numbers
    if evaluate_expression([eval(f"{numbers[0]}{ops[0]}{numbers[1]}"), numbers[2]], [ops[1]]) == 10:
        print(['(', numbers[0], ops[0], numbers[1], ')', ops[1], numbers[2]])
    
    # Evaluate with parentheses around last two numbers
    if evaluate_expression([numbers[0], eval(f"{numbers[1]}{ops[1]}{numbers[2]}")], [ops[0]]) == 10:
        print([numbers[0], ops[0], '(', numbers[1], ops[1], numbers[2], ')'])