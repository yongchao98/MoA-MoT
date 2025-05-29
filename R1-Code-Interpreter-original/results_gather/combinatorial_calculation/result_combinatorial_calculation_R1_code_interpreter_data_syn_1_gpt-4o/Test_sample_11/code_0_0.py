from itertools import product

# Given numbers
numbers = [9, 6, 4, 6]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all possible combinations of operations
for ops in product(operations, repeat=3):
    # Generate expressions with different parenthesis placements
    expressions = [
        f"({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} ({numbers[2]} {ops[2]} {numbers[3]})",
        f"(({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]}) {ops[2]} {numbers[3]}",
        f"({numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} {numbers[2]})) {ops[2]} {numbers[3]}",
        f"{numbers[0]} {ops[0]} (({numbers[1]} {ops[1]} {numbers[2]}) {ops[2]} {numbers[3]})",
        f"{numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} ({numbers[2]} {ops[2]} {numbers[3]}))"
    ]
    
    # Evaluate each expression
    for expr in expressions:
        if evaluate_expression(expr) == 24:
            print(expr)
            break