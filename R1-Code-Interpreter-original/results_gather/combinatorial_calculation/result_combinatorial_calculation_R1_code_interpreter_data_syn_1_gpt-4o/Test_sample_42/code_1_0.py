from itertools import product

# Define the numbers and target
numbers = [9, 6, 1, 7, 1, 3]
target = 100

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all combinations of operations
for ops in product(operations, repeat=5):
    # Generate the expression with the current combination of operations
    expr = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}{ops[4]}{numbers[5]}"
    
    # Evaluate the expression
    result = evaluate_expression(expr)
    
    # Check if the result matches the target
    if result == target:
        print(f"Expression: {expr} = {result}")
        break