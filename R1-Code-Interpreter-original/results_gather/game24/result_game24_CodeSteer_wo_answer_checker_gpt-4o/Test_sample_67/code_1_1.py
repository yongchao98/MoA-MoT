import itertools
import operator

# Define the numbers and target
numbers = [1, 2, 5, 7]
target = 24

# Define possible operations
operations = [operator.add, operator.sub, operator.mul, operator.truediv]
operation_symbols = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    """Evaluate the expression given numbers and operations."""
    try:
        result = nums[0]
        expression = str(nums[0])
        for i in range(3):
            result = ops[i](result, nums[i + 1])
            expression += f" {operation_symbols[operations.index(ops[i])]} {nums[i + 1]}"
        return result, expression
    except ZeroDivisionError:
        return None, None

def find_expression():
    """Find an expression that evaluates to the target."""
    for num_perm in itertools.permutations(numbers):
        for ops in itertools.product(operations, repeat=3):
            result, expression = evaluate_expression(num_perm, ops)
            if result == target:
                return expression
    return None

# Find and print the expression
expression = find_expression()
if expression:
    print(f"<<<{expression} = {target}>>>")
else:
    print("No solution found.")