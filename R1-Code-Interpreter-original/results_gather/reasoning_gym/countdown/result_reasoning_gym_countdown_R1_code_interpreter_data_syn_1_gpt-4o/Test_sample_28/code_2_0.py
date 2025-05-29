# Let's try to find a combination of operations that results in 259
from itertools import permutations, product

numbers = [83, 12, 68, 49, 47]
target = 259
operations = ['+', '-', '*', '/']

# Generate all permutations of numbers and operations
for num_perm in permutations(numbers):
    for ops in product(operations, repeat=4):
        # Construct the expression
        expr = f"({num_perm[0]} {ops[0]} ({num_perm[1]} {ops[1]} {num_perm[2]})) {ops[2]} {num_perm[3]} {ops[3]} {num_perm[4]}"
        try:
            # Evaluate the expression
            if eval(expr) == target:
                print(expr)
                break
        except ZeroDivisionError:
            continue