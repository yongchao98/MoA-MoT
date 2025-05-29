from itertools import product

# Define the numbers and target
numbers = [9, 7, 1, 1, 3]
target = 78
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(numbers, ops):
    expression = f"{numbers[0]}"
    for i in range(1, len(numbers)):
        expression += f"{ops[i-1]}{numbers[i]}"
    try:
        return eval(expression)
    except ZeroDivisionError:
        return None

# Try all combinations of operations
for ops in product(operations, repeat=4):
    if evaluate_expression(numbers, ops) == target:
        print([numbers[0], ops[0], numbers[1], ops[1], numbers[2], ops[2], numbers[3], ops[3], numbers[4]])
        break