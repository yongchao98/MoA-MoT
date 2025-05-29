from itertools import permutations, product

numbers = [51, 46, 36, 38, 21]
target = 537

# Try all permutations of the numbers
for perm in permutations(numbers):
    # Try all combinations of operations
    for ops in product('+-*/', repeat=4):
        # Construct the expression
        expr = f"{perm[0]}{ops[0]}{perm[1]}{ops[1]}{perm[2]}{ops[2]}{perm[3]}{ops[3]}{perm[4]}"
        try:
            # Evaluate the expression
            if eval(expr) == target:
                print(expr)
                break
        except ZeroDivisionError:
            continue