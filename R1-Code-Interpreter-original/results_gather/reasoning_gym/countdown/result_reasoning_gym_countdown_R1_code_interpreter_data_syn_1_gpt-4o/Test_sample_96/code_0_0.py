from itertools import permutations, product

numbers = [59, 74, 35, 58]
target = 583
operations = ['+', '-', '*', '/']

# Try all permutations of numbers and operations
for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        # Construct the expression
        expr1 = f"({num_perm[0]} {ops[0]} {num_perm[1]}) {ops[1]} {num_perm[2]} {ops[2]} {num_perm[3]}"
        expr2 = f"{num_perm[0]} {ops[0]} ({num_perm[1]} {ops[1]} {num_perm[2]}) {ops[2]} {num_perm[3]}"
        expr3 = f"{num_perm[0]} {ops[0]} {num_perm[1]} {ops[1]} ({num_perm[2]} {ops[2]} {num_perm[3]})"
        expr4 = f"({num_perm[0]} {ops[0]} {num_perm[1]} {ops[1]} {num_perm[2]}) {ops[2]} {num_perm[3]}"
        expr5 = f"{num_perm[0]} {ops[0]} ({num_perm[1]} {ops[1]} {num_perm[2]} {ops[2]} {num_perm[3]})"
        
        for expr in [expr1, expr2, expr3, expr4, expr5]:
            try:
                if eval(expr) == target:
                    print(expr)
            except ZeroDivisionError:
                continue