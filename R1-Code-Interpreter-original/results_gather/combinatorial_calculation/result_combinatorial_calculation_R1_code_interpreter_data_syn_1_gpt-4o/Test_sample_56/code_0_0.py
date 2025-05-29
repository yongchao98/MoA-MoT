from itertools import product

# Define the numbers and target
numbers = [9, 8, 8]
target = 10

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all possible combinations of operations
for ops in product(operations, repeat=2):
    # Generate expressions with different parenthesis placements
    expressions = [
        f"({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]}",
        f"{numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} {numbers[2]})",
        f"({numbers[0]} {ops[0]} {numbers[1]} {ops[1]} {numbers[2]})"
    ]
    
    # Evaluate each expression
    for expr in expressions:
        if evaluate_expression(expr) == target:
            print(expr)