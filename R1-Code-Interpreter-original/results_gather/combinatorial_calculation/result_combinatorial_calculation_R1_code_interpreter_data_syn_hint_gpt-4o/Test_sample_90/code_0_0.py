from itertools import product

# Define the numbers and target
numbers = [7, 3, 2, 9, 8]
target = 78
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all possible combinations of operations
for ops in product(operations, repeat=4):
    # Create the expression with the current combination of operations
    expr = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}"
    
    # Evaluate the expression
    if evaluate_expression(expr) == target:
        print(expr)
        break