from sympy import symbols, Eq, solve, sympify
from itertools import product

# Define the numbers and operations
numbers = [6, 2, 5, 2, 1, 4]
operations = ['+', '-', '*', '/']

# Function to evaluate an expression using sympy
def evaluate_expression(expr):
    try:
        return sympify(expr).evalf()
    except (ZeroDivisionError, ValueError):
        return None

# Generate all possible combinations of operations
for ops in product(operations, repeat=5):
    # Try different parenthesis placements
    expressions = [
        f"(({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]}) {ops[2]} ({numbers[3]} {ops[3]} ({numbers[4]} {ops[4]} {numbers[5]}))",
        f"({numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} {numbers[2]})) {ops[2]} ({numbers[3]} {ops[3]} ({numbers[4]} {ops[4]} {numbers[5]}))",
        f"({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} (({numbers[2]} {ops[2]} {numbers[3]}) {ops[3]} ({numbers[4]} {ops[4]} {numbers[5]}))",
        f"({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} ({numbers[2]} {ops[2]} ({numbers[3]} {ops[3]} ({numbers[4]} {ops[4]} {numbers[5]})))",
        f"(({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]}) {ops[2]} {numbers[3]} {ops[3]} ({numbers[4]} {ops[4]} {numbers[5]})",
    ]
    
    for expr in expressions:
        if evaluate_expression(expr) == 100:
            print(expr)
            break