import itertools

numbers = [28, 61, 63, 14, 75, 2]
target = 117

# Define a function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all permutations of the numbers
for perm in itertools.permutations(numbers):
    # Generate all possible combinations of operations
    for ops in itertools.product('+-*/', repeat=5):
        # Construct the expression
        expr = f"{perm[0]}{ops[0]}{perm[1]}{ops[1]}{perm[2]}{ops[2]}{perm[3]}{ops[3]}{perm[4]}{ops[4]}{perm[5]}"
        # Evaluate the expression
        if evaluate_expression(expr) == target:
            print(expr)
            break