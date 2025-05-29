from itertools import permutations
import operator

# Define the numbers and target
numbers = [1, 4, 6, 12]
target = 24

# Define possible operations
ops = [operator.add, operator.sub, operator.mul, operator.truediv]
op_symbols = ['+', '-', '*', '/']

# Function to evaluate expression
def evaluate_expression(nums, ops):
    try:
        result = ops[0](nums[0], nums[1])
        result = ops[1](result, nums[2])
        result = ops[2](result, nums[3])
        return result
    except ZeroDivisionError:
        return None

# Try all permutations of numbers and operations
for num_perm in permutations(numbers):
    for op_perm in permutations(ops, 3):
        if evaluate_expression(num_perm, op_perm) == target:
            # Construct the expression string
            expr = f"(({num_perm[0]} {op_symbols[ops.index(op_perm[0])]} {num_perm[1]}) {op_symbols[ops.index(op_perm[1])]} {num_perm[2]}) {op_symbols[ops.index(op_perm[2])]} {num_perm[3]}"
            print(expr)
            break