from itertools import product

# Given numbers
numbers = [9, 1, 7, 9, 7]
# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all combinations of operations and parentheses
for ops in product(operations, repeat=4):
    # Generate possible expressions with parentheses
    expressions = [
        f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}{numbers[3]}){ops[3]}{numbers[4]}",
        f"({numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}){ops[2]}({numbers[3]}{ops[3]}{numbers[4]})",
        f"{numbers[0]}{ops[0]}({numbers[1]}{ops[1]}({numbers[2]}{ops[2]}{numbers[3]})){ops[3]}{numbers[4]}",
        f"{numbers[0]}{ops[0]}({numbers[1]}{ops[1]}{numbers[2]}){ops[2]}({numbers[3]}{ops[3]}{numbers[4]})",
        f"({numbers[0]}{ops[0]}({numbers[1]}{ops[1]}{numbers[2]})){ops[2]}({numbers[3]}{ops[3]}{numbers[4]})"
    ]
    
    for expr in expressions:
        if evaluate_expression(expr) == 78:
            print(expr)
            break