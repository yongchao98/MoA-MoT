from itertools import product

# Define the numbers and target
numbers = [3, 6, 8, 2]
target = 24

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression
def evaluate_expression(nums, ops):
    try:
        # Construct the expression
        expression = f"({nums[0]} {ops[0]} {nums[1]}) {ops[1]} ({nums[2]} {ops[2]} {nums[3]})"
        # Evaluate the expression
        return eval(expression)
    except ZeroDivisionError:
        return None

# Try all combinations of operations
for ops in product(operations, repeat=3):
    result = evaluate_expression(numbers, ops)
    if result == target:
        print(f"Solution found: ({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} ({numbers[2]} {ops[2]} {numbers[3]})")
        break