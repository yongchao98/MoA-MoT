from itertools import product

# The numbers in the given order
numbers = [7, 4, 5, 5, 5, 5]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression safely
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all combinations of operations
for ops in product(operations, repeat=5):
    # Generate the expression without parentheses
    expr = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}{ops[4]}{numbers[5]}"
    if evaluate_expression(expr) == 100:
        print(f"Solution without parentheses: {expr}")
    
    # Try adding parentheses in different places
    for i in range(1, 6):
        for j in range(i+1, 7):
            expr_with_parens = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}({numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}){ops[4]}{numbers[5]}"
            if evaluate_expression(expr_with_parens) == 100:
                print(f"Solution with parentheses: {expr_with_parens}")