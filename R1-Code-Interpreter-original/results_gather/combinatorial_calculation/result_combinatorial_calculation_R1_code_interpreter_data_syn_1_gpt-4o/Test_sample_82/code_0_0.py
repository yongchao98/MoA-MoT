from itertools import product

# Define the numbers and target
numbers = [3, 1, 9, 6]
target = 24

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all combinations of operations and parentheses
for ops in product(operations, repeat=3):
    # Generate possible expressions with parentheses
    expressions = [
        f"({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} ({numbers[2]} {ops[2]} {numbers[3]})",
        f"(({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]}) {ops[2]} {numbers[3]}",
        f"({numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} {numbers[2]})) {ops[2]} {numbers[3]}",
        f"{numbers[0]} {ops[0]} (({numbers[1]} {ops[1]} {numbers[2]}) {ops[2]} {numbers[3]})",
        f"{numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} ({numbers[2]} {ops[2]} {numbers[3]}))"
    ]
    
    # Evaluate each expression
    for expr in expressions:
        if evaluate_expression(expr) == target:
            print(expr)
            break