import itertools

# Define the numbers and the target
numbers = [63, 26, 52, 60, 28]
target = 889

# Define the operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all permutations of numbers and operations
for num_perm in itertools.permutations(numbers):
    for ops in itertools.product(operations, repeat=4):
        # Create the expression
        expr = f"(({num_perm[0]} {ops[0]} {num_perm[1]}) {ops[1]} ({num_perm[2]} {ops[2]} {num_perm[3]})) {ops[3]} {num_perm[4]}"
        if evaluate_expression(expr) == target:
            print(expr)
            break