from itertools import permutations, product

# Define the numbers and target
numbers = [8, 5, 8, 8]
target = 24
operations = ['+', '-', '*', '/']

# Function to evaluate an expression safely
def safe_eval(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all permutations of numbers
number_permutations = permutations(numbers)

# Generate all combinations of operations
operation_combinations = product(operations, repeat=3)

# Try all combinations of numbers, operations, and parentheses
for num_perm in number_permutations:
    for ops in operation_combinations:
        # Generate all possible expressions with parentheses
        expressions = [
            f"({num_perm[0]} {ops[0]} {num_perm[1]}) {ops[1]} ({num_perm[2]} {ops[2]} {num_perm[3]})",
            f"(({num_perm[0]} {ops[0]} {num_perm[1]}) {ops[1]} {num_perm[2]}) {ops[2]} {num_perm[3]}",
            f"({num_perm[0]} {ops[0]} ({num_perm[1]} {ops[1]} {num_perm[2]})) {ops[2]} {num_perm[3]}",
            f"{num_perm[0]} {ops[0]} (({num_perm[1]} {ops[1]} {num_perm[2]}) {ops[2]} {num_perm[3]})",
            f"{num_perm[0]} {ops[0]} ({num_perm[1]} {ops[1]} ({num_perm[2]} {ops[2]} {num_perm[3]}))"
        ]
        
        # Evaluate each expression
        for expr in expressions:
            if safe_eval(expr) == target:
                print(expr)
                break