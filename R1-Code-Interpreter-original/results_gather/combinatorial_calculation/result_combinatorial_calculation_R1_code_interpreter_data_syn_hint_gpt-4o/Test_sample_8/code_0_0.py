from itertools import product

# Define the numbers and target
numbers = [9, 7, 6]
target = 10

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression
def evaluate_expression(nums, ops):
    try:
        # Construct the expression
        expression = f"{nums[0]} {ops[0]} {nums[1]} {ops[1]} {nums[2]}"
        # Evaluate the expression
        return eval(expression)
    except ZeroDivisionError:
        return None

# Try all combinations of operations
for ops in product(operations, repeat=2):
    if evaluate_expression(numbers, ops) == target:
        print([numbers[0], ops[0], numbers[1], ops[1], numbers[2]])
        break