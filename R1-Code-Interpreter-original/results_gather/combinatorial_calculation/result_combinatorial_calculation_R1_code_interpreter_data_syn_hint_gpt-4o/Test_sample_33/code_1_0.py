from itertools import product

# Given numbers
numbers = [7, 8, 2, 4, 4]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression with given operations
def evaluate_expression(numbers, ops):
    expr = f"{numbers[0]}"
    for i in range(1, len(numbers)):
        expr += f"{ops[i-1]}{numbers[i]}"
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all combinations of operations
for ops in product(operations, repeat=4):
    if evaluate_expression(numbers, ops) == 69:
        print([numbers[0], ops[0], numbers[1], ops[1], numbers[2], ops[2], numbers[3], ops[3], numbers[4]])
        break