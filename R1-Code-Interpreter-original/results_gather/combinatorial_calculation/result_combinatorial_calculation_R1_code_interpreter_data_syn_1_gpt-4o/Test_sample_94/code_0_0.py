from itertools import product

# The numbers in the given order
numbers = [8, 3, 7, 1, 9]

# The target result
target = 78

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all combinations of operations
for ops in product(operations, repeat=4):
    # Generate the expression without parentheses
    expr = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}"
    if evaluate_expression(expr) == target:
        print(f"Solution without parentheses: {expr}")
    
    # Try adding parentheses in different places
    for i in range(4):
        for j in range(i+1, 5):
            expr_with_parens = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}({numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]})"
            if evaluate_expression(expr_with_parens) == target:
                print(f"Solution with parentheses: {expr_with_parens}")