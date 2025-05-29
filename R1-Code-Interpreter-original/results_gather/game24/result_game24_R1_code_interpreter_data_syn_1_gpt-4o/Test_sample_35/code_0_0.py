import itertools

# Define the numbers and target
numbers = [1, 2, 6, 12]
target = 24

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all permutations of numbers
number_permutations = itertools.permutations(numbers)

# Generate all combinations of operations
operation_combinations = itertools.product(operations, repeat=3)

# Try all combinations of numbers and operations
for num_perm in number_permutations:
    for ops in operation_combinations:
        # Generate possible expressions with different parenthesis placements
        expressions = [
            f"(({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]}",
            f"({num_perm[0]}{ops[0]}({num_perm[1]}{ops[1]}{num_perm[2]})){ops[2]}{num_perm[3]}",
            f"{num_perm[0]}{ops[0]}(({num_perm[1]}{ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]})",
            f"{num_perm[0]}{ops[0]}({num_perm[1]}{ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]}))",
            f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})"
        ]
        
        # Evaluate each expression
        for expr in expressions:
            if evaluate_expression(expr) == target:
                print(expr)
                break