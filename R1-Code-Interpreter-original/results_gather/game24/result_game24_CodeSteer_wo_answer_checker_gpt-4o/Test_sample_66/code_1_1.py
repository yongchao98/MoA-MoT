import itertools
import operator

# Define the numbers and target
numbers = [5, 7, 10, 12]
target = 24

# Define possible operations
operations = [operator.add, operator.sub, operator.mul, operator.truediv]
operation_symbols = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    """Evaluate the expression with given numbers and operations."""
    try:
        # Apply operations in the order of precedence
        result = ops[0](nums[0], nums[1])
        result = ops[1](result, nums[2])
        result = ops[2](result, nums[3])
        return result
    except ZeroDivisionError:
        return None

def find_expression():
    """Find a valid expression that evaluates to the target."""
    for num_perm in itertools.permutations(numbers):
        for ops in itertools.product(operations, repeat=3):
            if evaluate_expression(num_perm, ops) == target:
                # Format the expression with the correct symbols
                expr = f"(({num_perm[0]}{operation_symbols[operations.index(ops[0])]}{num_perm[1]})" \
                       f"{operation_symbols[operations.index(ops[1])]}{num_perm[2]})" \
                       f"{operation_symbols[operations.index(ops[2])]}{num_perm[3]}"
                return expr
    return None

# Find and print the expression
expression = find_expression()
if expression:
    print(f"<<<{expression}>>>")
else:
    print("No valid expression found.")