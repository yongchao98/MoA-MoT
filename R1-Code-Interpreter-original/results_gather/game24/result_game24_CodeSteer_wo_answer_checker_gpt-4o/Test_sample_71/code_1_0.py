import itertools
import operator

# Define the numbers and target
numbers = [8, 9, 11, 12]
target = 24

# Define possible operations
operations = [operator.add, operator.sub, operator.mul, operator.truediv]
operation_symbols = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    """Evaluate the expression with given numbers and operations."""
    try:
        result = ops[0](nums[0], ops[1](nums[1], ops[2](nums[2], nums[3])))
        return result
    except ZeroDivisionError:
        return None

def find_expression():
    """Find an expression that evaluates to the target."""
    for num_perm in itertools.permutations(numbers):
        for ops in itertools.product(operations, repeat=3):
            result = evaluate_expression(num_perm, ops)
            if result == target:
                # Construct the expression string
                expr = f"(({num_perm[0]}{operation_symbols[operations.index(ops[0])]}({num_perm[1]}{operation_symbols[operations.index(ops[1])]}({num_perm[2]}{operation_symbols[operations.index(ops[2])]}{num_perm[3]}))))"
                return expr
    return None

# Find and print the expression
expression = find_expression()
if expression:
    print(f"<<<{expression} = {target}>>>")
else:
    print("No valid expression found.")