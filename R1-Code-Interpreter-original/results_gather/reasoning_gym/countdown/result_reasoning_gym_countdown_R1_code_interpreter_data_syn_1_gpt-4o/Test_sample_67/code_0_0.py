import itertools

# Define the numbers and target
numbers = [96, 65, 77, 23, 57]
target = 204

# Define the operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all permutations of numbers
for num_perm in itertools.permutations(numbers):
    # Generate all combinations of operations
    for ops in itertools.product(operations, repeat=4):
        # Create the expression
        expr = f"({num_perm[0]} {ops[0]} {num_perm[1]}) {ops[1]} ({num_perm[2]} {ops[2]} {num_perm[3]}) {ops[3]} {num_perm[4]}"
        # Evaluate the expression
        if evaluate_expression(expr) == target:
            print(expr)
            break