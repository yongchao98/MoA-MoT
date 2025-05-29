import itertools

# Define the numbers and target
numbers = [1, 2, 4, 12]
target = 24
operations = ['+', '-', '*', '/']

# Function to evaluate an expression safely
def safe_eval(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all permutations of the numbers
for num_perm in itertools.permutations(numbers):
    # Generate all combinations of operations
    for ops in itertools.product(operations, repeat=3):
        # Generate all possible expressions with different groupings
        expressions = [
            f"(({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]}",
            f"({num_perm[0]}{ops[0]}({num_perm[1]}{ops[1]}{num_perm[2]})){ops[2]}{num_perm[3]}",
            f"{num_perm[0]}{ops[0]}(({num_perm[1]}{ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]})",
            f"{num_perm[0]}{ops[0]}({num_perm[1]}{ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]}))",
            f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})"
        ]
        
        # Evaluate each expression
        for expr in expressions:
            if safe_eval(expr) == target:
                print(f"<<<{expr}={target}>>>")
                break