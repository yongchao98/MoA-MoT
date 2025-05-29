from itertools import product

# Define the numbers and target
numbers = [4, 4, 2]
target = 10

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all combinations of operations and parentheses
for ops in product(operations, repeat=2):
    expr1 = f"({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]}"
    expr2 = f"{numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} {numbers[2]})"
    expr3 = f"({numbers[0]} {ops[0]} {numbers[1]} {ops[1]} {numbers[2]})"
    
    for expr in [expr1, expr2, expr3]:
        if evaluate_expression(expr) == target:
            print(expr)
            break