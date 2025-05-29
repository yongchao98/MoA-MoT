from itertools import product

# Given numbers
numbers = [6, 9, 6, 7, 2, 8]
# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression with parentheses
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Recursive function to generate expressions
def generate_expressions(numbers, ops, index=0, current_expr=""):
    if index == len(numbers) - 1:
        # Base case: evaluate the final expression
        if evaluate_expression(current_expr) == 100:
            return current_expr
        return None
    
    # Recursive case: try each operation
    for op in ops:
        new_expr = f"({current_expr}{op}{numbers[index+1]})"
        result = generate_expressions(numbers, ops, index+1, new_expr)
        if result:
            return result
    return None

# Start with the first number
initial_expr = str(numbers[0])
# Find a valid expression
valid_expression = generate_expressions(numbers, operations, 0, initial_expr)

if valid_expression:
    print(valid_expression)
else:
    print("No valid expression found.")