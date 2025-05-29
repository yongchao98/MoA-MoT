from itertools import product

# Define the numbers and target
numbers = [6, 2, 5, 2, 1, 4]
target = 100

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression safely
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all possible combinations of operations
for ops in product(operations, repeat=5):
    # Create the expression string
    expr = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}{ops[4]}{numbers[5]}"
    
    # Evaluate the expression
    if evaluate_expression(expr) == target:
        # If the expression evaluates to the target, print the expression
        print(expr)
        break