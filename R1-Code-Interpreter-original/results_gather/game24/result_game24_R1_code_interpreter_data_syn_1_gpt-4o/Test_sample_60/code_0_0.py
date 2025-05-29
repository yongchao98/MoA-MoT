import itertools

# Define the numbers and target
numbers = [2, 8, 8, 13]
target = 24

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
    for ops in itertools.product(operations, repeat=3):
        # Generate all possible parenthesized expressions
        expressions = [
            f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})",
            f"(({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]}",
            f"({num_perm[0]}{ops[0]}({num_perm[1]}{ops[1]}{num_perm[2]})){ops[2]}{num_perm[3]}",
            f"{num_perm[0]}{ops[0]}(({num_perm[1]}{ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]})",
            f"{num_perm[0]}{ops[0]}({num_perm[1]}{ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]}))"
        ]
        # Evaluate each expression
        for expr in expressions:
            if evaluate_expression(expr) == target:
                print(expr)
                break