from itertools import product

# Define the numbers and target
numbers = [7, 6, 1, 2]
target = 24

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(nums, ops):
    try:
        # Construct the expression
        expression = f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})"
        # Evaluate the expression
        result = eval(expression)
        return result
    except ZeroDivisionError:
        return None

# Try all combinations of operations
for ops in product(operations, repeat=3):
    if evaluate_expression(numbers, ops) == target:
        # If a valid expression is found, print it
        expression = f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}{numbers[3]})"
        print(expression)
        break