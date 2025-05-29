import itertools
import operator

# Define the numbers and target
numbers = [3, 8, 11, 11]
target = 24

# Define possible operations
operations = [operator.add, operator.sub, operator.mul, operator.truediv]

# Function to evaluate an expression
def evaluate_expression(nums, ops):
    try:
        result = ops[0](nums[0], nums[1])
        result = ops[1](result, nums[2])
        result = ops[2](result, nums[3])
        return result
    except ZeroDivisionError:
        return None

# Try all permutations of numbers and combinations of operations
for num_perm in itertools.permutations(numbers):
    for ops in itertools.product(operations, repeat=3):
        if evaluate_expression(num_perm, ops) == target:
            expression = f"(({num_perm[0]} {ops[0].__name__} {num_perm[1]}) {ops[1].__name__} {num_perm[2]}) {ops[2].__name__} {num_perm[3]}"
            print(expression)
            break