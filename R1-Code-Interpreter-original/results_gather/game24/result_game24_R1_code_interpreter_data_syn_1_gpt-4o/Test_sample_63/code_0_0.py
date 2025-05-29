from itertools import permutations, product

# Define the numbers and target
numbers = [2, 2, 8, 12]
target = 24

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(nums, ops):
    try:
        # Construct the expression
        expression = f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]}))"
        # Evaluate the expression
        result = eval(expression)
        # Check if the result is equal to the target
        if result == target:
            return expression
    except ZeroDivisionError:
        return None
    return None

# Try all permutations of numbers and operations
for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        expr = evaluate_expression(num_perm, ops)
        if expr:
            print(expr)
            break