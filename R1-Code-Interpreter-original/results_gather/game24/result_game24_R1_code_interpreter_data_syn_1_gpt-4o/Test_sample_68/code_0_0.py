# Define the numbers
numbers = [6, 8, 10, 12]

# Try to find a combination that results in 24
# We will use a brute force approach to try different combinations

from itertools import permutations, product

# Define the operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all permutations of numbers and operations
for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        # Form expressions with different parenthesis placements
        expressions = [
            f"(({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]}",
            f"({num_perm[0]}{ops[0]}({num_perm[1]}{ops[1]}{num_perm[2]})){ops[2]}{num_perm[3]}",
            f"{num_perm[0]}{ops[0]}(({num_perm[1]}{ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]})",
            f"{num_perm[0]}{ops[0]}({num_perm[1]}{ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]}))",
            f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})"
        ]
        
        # Evaluate each expression
        for expr in expressions:
            if evaluate_expression(expr) == 24:
                print(expr)
                break