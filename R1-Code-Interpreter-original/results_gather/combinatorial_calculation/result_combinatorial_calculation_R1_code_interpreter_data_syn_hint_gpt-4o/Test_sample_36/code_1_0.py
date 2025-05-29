from itertools import product, permutations

# Given numbers
numbers = [9, 6, 5, 7, 7]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all possible expressions
def generate_expressions(numbers, operations):
    for ops in product(operations, repeat=4):
        # Generate expression without parentheses
        expr = f"{numbers[0]} {ops[0]} {numbers[1]} {ops[1]} {numbers[2]} {ops[2]} {numbers[3]} {ops[3]} {numbers[4]}"
        yield expr
        # Generate expressions with parentheses
        for i in range(4):
            for j in range(i+2, 6):
                expr_with_parens = expr[:i*4] + '(' + expr[i*4:j*4-1] + ')' + expr[j*4-1:]
                yield expr_with_parens

# Find the correct expression
for expr in generate_expressions(numbers, operations):
    if evaluate_expression(expr) == 69:
        print(expr)
        break