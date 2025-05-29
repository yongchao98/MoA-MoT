from itertools import permutations, product

def find_combination():
    numbers = [3, 3, 61, 36, 75]
    target = 654
    operations = ['+', '-', '*', '/']
    
    # Try all permutations of numbers
    for num_perm in permutations(numbers):
        # Try all combinations of operations
        for ops in product(operations, repeat=4):
            # Construct the expression
            expr = f"{num_perm[0]} {ops[0]} {num_perm[1]} {ops[1]} {num_perm[2]} {ops[2]} {num_perm[3]} {ops[3]} {num_perm[4]}"
            try:
                # Evaluate the expression
                if eval(expr) == target:
                    return expr
            except ZeroDivisionError:
                continue
    
    return "No combination found"

print(find_combination())