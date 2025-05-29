from itertools import permutations
import operator

# Define the numbers and target
numbers = [7, 9, 11, 11]
target = 24

# Define possible operations
operations = [operator.add, operator.sub, operator.mul, operator.truediv]
operation_symbols = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    """Evaluate the expression given numbers and operations."""
    result = nums[0]
    expression = str(nums[0])
    for i in range(3):
        result = ops[i](result, nums[i + 1])
        expression = f"({expression} {operation_symbols[operations.index(ops[i])]} {nums[i + 1]})"
    return result, expression

def find_solution():
    """Find a solution that results in the target value."""
    for perm in permutations(numbers):
        for ops in permutations(operations, 3):
            try:
                result, expression = evaluate_expression(perm, ops)
                if abs(result - target) < 1e-6:  # Account for floating-point precision
                    return expression
            except ZeroDivisionError:
                continue
    return None

solution = find_solution()
if solution:
    print(f"<<<{solution} = {target}>>>")
else:
    print("No solution found.")