import itertools

# Define the numbers and target
numbers = [2, 6, 6, 8]
target = 24

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(nums, ops):
    try:
        # Construct the expression
        expression = f"(({nums[0]} {ops[0]} {nums[1]}) {ops[1]} ({nums[2]} {ops[2]} {nums[3]}))"
        # Evaluate the expression
        result = eval(expression)
        return result
    except ZeroDivisionError:
        return None

# Iterate over all permutations of numbers and operations
for num_perm in itertools.permutations(numbers):
    for ops in itertools.product(operations, repeat=3):
        result = evaluate_expression(num_perm, ops)
        if result == target:
            # Construct the expression string
            expression = f"(({num_perm[0]} {ops[0]} {num_perm[1]}) {ops[1]} ({num_perm[2]} {ops[2]} {num_perm[3]}))"
            print(expression)
            break