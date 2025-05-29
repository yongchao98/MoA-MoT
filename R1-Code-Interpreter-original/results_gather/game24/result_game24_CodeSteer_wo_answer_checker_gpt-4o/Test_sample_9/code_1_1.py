from itertools import permutations
from operator import add, sub, mul, truediv

# Define the numbers and target
numbers = [2, 2, 3, 8]
target = 24

# Define the operations
operations = [add, sub, mul, truediv]
operation_symbols = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    """Evaluate the expression given numbers and operations."""
    result = nums[0]
    expression = str(nums[0])
    for i in range(3):
        result = ops[i](result, nums[i + 1])
        expression = f"({expression} {operation_symbols[ops[i]]} {nums[i + 1]})"
    return result, expression

def find_expression():
    """Find an expression that evaluates to the target."""
    for num_perm in permutations(numbers):
        for ops in permutations(operations, 3):
            try:
                result, expression = evaluate_expression(num_perm, ops)
                if abs(result - target) < 1e-6:  # Check for floating-point precision
                    return expression
            except ZeroDivisionError:
                continue
    return None

# Find and print the expression
expression = find_expression()
if expression:
    print(f"<<<{expression}>>>")
else:
    print("No valid expression found.")