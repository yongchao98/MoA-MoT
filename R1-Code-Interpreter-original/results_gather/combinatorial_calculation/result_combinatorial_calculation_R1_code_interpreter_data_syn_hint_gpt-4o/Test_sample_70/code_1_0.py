from itertools import permutations, product

# Define the numbers and target
numbers = [6, 5, 2, 2]
target = 24

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression
def evaluate_expression(nums, ops):
    try:
        # Construct the expression
        expression = f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})"
        # Evaluate the expression
        return eval(expression)
    except ZeroDivisionError:
        return None

# Try all permutations of numbers and combinations of operations
for nums in permutations(numbers):
    for ops in product(operations, repeat=3):
        if evaluate_expression(nums, ops) == target:
            result = [nums[0], ops[0], nums[1], ops[1], '(', nums[2], ops[2], nums[3], ')']
            print(result)
            break