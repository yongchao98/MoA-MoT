from itertools import permutations

numbers = [51, 4, 32, 13]
target = 50

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all permutations of the numbers
for perm in permutations(numbers):
    # Try different operations between the numbers
    for ops in permutations(['+', '-', '*', '/'], 3):
        # Form the expression
        expr1 = f"({perm[0]} {ops[0]} {perm[1]}) {ops[1]} {perm[2]} {ops[2]} {perm[3]}"
        expr2 = f"{perm[0]} {ops[0]} ({perm[1]} {ops[1]} {perm[2]}) {ops[2]} {perm[3]}"
        expr3 = f"{perm[0]} {ops[0]} {perm[1]} {ops[1]} ({perm[2]} {ops[2]} {perm[3]})"
        expr4 = f"({perm[0]} {ops[0]} {perm[1]} {ops[1]} {perm[2]}) {ops[2]} {perm[3]}"
        expr5 = f"{perm[0]} {ops[0]} ({perm[1]} {ops[1]} ({perm[2]} {ops[2]} {perm[3]}))"
        
        # Evaluate the expressions
        for expr in [expr1, expr2, expr3, expr4, expr5]:
            if evaluate_expression(expr) == target:
                print(expr)