from sympy import symbols, Eq, solve
from itertools import product

# Define the numbers and operations
numbers = [9, 6, 5, 7, 7]
operations = ['+', '-', '*', '/']

# Define symbolic variables
a, b, c, d, e = symbols('a b c d e')

# Function to evaluate expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all possible expressions
def generate_expressions(numbers, operations):
    for ops in product(operations, repeat=4):
        # Create expression with symbolic variables
        expr = f"{a} {ops[0]} {b} {ops[1]} {c} {ops[2]} {d} {ops[3]} {e}"
        # Substitute numbers into the expression
        expr_sub = expr.replace('a', str(numbers[0])).replace('b', str(numbers[1])).replace('c', str(numbers[2])).replace('d', str(numbers[3])).replace('e', str(numbers[4]))
        yield expr_sub

# Find the correct expression
for expr in generate_expressions(numbers, operations):
    if evaluate_expression(expr) == 69:
        print(expr)
        break