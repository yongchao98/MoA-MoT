from itertools import permutations
import operator

# Define the numbers and target
numbers = [1, 3, 5, 12]
target = 24

# Define possible operations
operations = [operator.add, operator.sub, operator.mul, operator.truediv]
operation_symbols = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    """Evaluate the expression with given numbers and operations."""
    try:
        result = ops[0](nums[0], nums[1])
        result = ops[1](result, nums[2])
        result = ops[2](result, nums[3])
        return result
    except ZeroDivisionError:
        return None

def find_solution():
    """Find a solution that results in the target value."""
    for num_perm in permutations(numbers):
        for ops in permutations(operations, 3):
            if evaluate_expression(num_perm, ops) == target:
                # Construct the expression string
                expr = f"(({num_perm[0]} {operation_symbols[operations.index(ops[0])]} {num_perm[1]}) {operation_symbols[operations.index(ops[1])]} {num_perm[2]}) {operation_symbols[operations.index(ops[2])]} {num_perm[3]}"
                return expr
    return None

# Find and print the solution
solution = find_solution()
if solution:
    print(f"<<<{solution} = {target}>>>")
else:
    print("No solution found.")