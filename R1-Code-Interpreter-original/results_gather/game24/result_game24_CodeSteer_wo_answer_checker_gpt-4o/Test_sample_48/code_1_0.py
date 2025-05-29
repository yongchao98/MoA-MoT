from itertools import permutations, product
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
        try:
            result = ops[i](result, nums[i + 1])
            expression = f"({expression} {operation_symbols[operations.index(ops[i])]} {nums[i + 1]})"
        except ZeroDivisionError:
            return None, None
    return result, expression

def find_solution():
    """Find a solution that results in the target value."""
    for perm in permutations(numbers):
        for ops in product(operations, repeat=3):
            result, expression = evaluate_expression(perm, ops)
            if result is not None and abs(result - target) < 1e-6:  # Account for floating-point precision
                return expression
    return None

solution = find_solution()
if solution:
    print(f"<<<{solution} = {target}>>>")
else:
    print("No solution found.")