from itertools import product

# Define the numbers and target
numbers = [7, 3, 2, 9, 8]
target = 78

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression safely
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all possible combinations of operations
for ops in product(operations, repeat=4):
    # Generate the expression without parentheses
    expr = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}"
    if evaluate_expression(expr) == target:
        print(f"Solution without parentheses: {expr}")
    
    # Generate expressions with parentheses
    # (a op b) op (c op d op e)
    expr = f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]})"
    if evaluate_expression(expr) == target:
        print(f"Solution with parentheses: {expr}")
    
    # a op (b op (c op (d op e)))
    expr = f"{numbers[0]}{ops[0]}({numbers[1]}{ops[1]}({numbers[2]}{ops[2]}({numbers[3]}{ops[3]}{numbers[4]})))"
    if evaluate_expression(expr) == target:
        print(f"Solution with nested parentheses: {expr}")