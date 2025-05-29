from itertools import product

# Define the numbers and target
numbers = [7, 4, 5, 7, 5, 7]
target = 100

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression safely
def safe_eval(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all combinations of operations
for ops in product(operations, repeat=5):
    # Generate the expression without parentheses
    expr = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}{ops[4]}{numbers[5]}"
    if safe_eval(expr) == target:
        print(f"Solution without parentheses: {expr}")
    
    # Try adding parentheses in different places
    for i in range(5):
        for j in range(i+2, 7):
            expr_with_parens = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}{ops[4]}{numbers[5]}"
            expr_with_parens = expr_with_parens[:i] + '(' + expr_with_parens[i:j] + ')' + expr_with_parens[j:]
            if safe_eval(expr_with_parens) == target:
                print(f"Solution with parentheses: {expr_with_parens}")