from itertools import product
from sympy import symbols, Eq, solve

# Define the numbers and operations
numbers = [4, 2, 6, 5]
operations = ['+', '-', '*', '/']

# Generate all possible combinations of operations
operation_combinations = list(product(operations, repeat=3))

# Define the symbols for the numbers
a, b, c, d = symbols('a b c d')

# Function to evaluate expressions
def evaluate_expression(expr):
    try:
        # Solve the equation to see if it equals 24
        return solve(Eq(expr, 24))
    except:
        return None

# Iterate over all combinations of operations
for ops in operation_combinations:
    # Generate expressions with different parenthesis placements
    expressions = [
        f"(({a} {ops[0]} {b}) {ops[1]} {c}) {ops[2]} {d}",
        f"({a} {ops[0]} ({b} {ops[1]} {c})) {ops[2]} {d}",
        f"{a} {ops[0]} (({b} {ops[1]} {c}) {ops[2]} {d})",
        f"{a} {ops[0]} ({b} {ops[1]} ({c} {ops[2]} {d}))",
        f"({a} {ops[0]} {b}) {ops[1]} ({c} {ops[2]} {d})"
    ]
    
    # Evaluate each expression
    for expr in expressions:
        expr_eval = evaluate_expression(eval(expr))
        if expr_eval:
            print(f"Expression: {expr} = 24")
            break