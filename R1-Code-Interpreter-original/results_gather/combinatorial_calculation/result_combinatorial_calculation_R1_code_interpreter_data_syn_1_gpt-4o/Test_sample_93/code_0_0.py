from itertools import product

# Given numbers
numbers = [5, 3, 9, 9, 6]
# Target result
target = 78
# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression with parentheses
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all possible combinations of operations
for ops in product(operations, repeat=4):
    # Generate all possible parenthesized expressions
    expressions = [
        f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}{numbers[3]}){ops[3]}{numbers[4]}",
        f"({numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}){ops[2]}({numbers[3]}{ops[3]}{numbers[4]})",
        f"({numbers[0]}{ops[0]}({numbers[1]}{ops[1]}{numbers[2]})){ops[2]}({numbers[3]}{ops[3]}{numbers[4]})",
        f"{numbers[0]}{ops[0]}(({numbers[1]}{ops[1]}{numbers[2]}){ops[2]}{numbers[3]}){ops[3]}{numbers[4]}",
        f"{numbers[0]}{ops[0]}({numbers[1]}{ops[1]}({numbers[2]}{ops[2]}{numbers[3]})){ops[3]}{numbers[4]}"
    ]
    
    # Evaluate each expression
    for expr in expressions:
        if evaluate_expression(expr) == target:
            print(expr)
            break