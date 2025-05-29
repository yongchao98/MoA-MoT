import itertools

# Define the numbers and target
numbers = [1, 8, 10, 13]
target = 24
operations = ['+', '-', '*', '/']

# Function to evaluate an expression safely
def safe_eval(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all permutations of numbers
for num_perm in itertools.permutations(numbers):
    # Generate all combinations of operations
    for ops in itertools.product(operations, repeat=3):
        # Test different groupings with parentheses
        expressions = [
            f"(({num_perm[0]} {ops[0]} {num_perm[1]}) {ops[1]} {num_perm[2]}) {ops[2]} {num_perm[3]}",
            f"({num_perm[0]} {ops[0]} ({num_perm[1]} {ops[1]} {num_perm[2]})) {ops[2]} {num_perm[3]}",
            f"{num_perm[0]} {ops[0]} (({num_perm[1]} {ops[1]} {num_perm[2]}) {ops[2]} {num_perm[3]})",
            f"{num_perm[0]} {ops[0]} ({num_perm[1]} {ops[1]} ({num_perm[2]} {ops[2]} {num_perm[3]}))",
            f"(({num_perm[0]} {ops[0]} {num_perm[1]}) {ops[1]} ({num_perm[2]} {ops[2]} {num_perm[3]}))"
        ]
        
        # Evaluate each expression
        for expr in expressions:
            if safe_eval(expr) == target:
                print(f"<<<{expr}>>>")
                exit()