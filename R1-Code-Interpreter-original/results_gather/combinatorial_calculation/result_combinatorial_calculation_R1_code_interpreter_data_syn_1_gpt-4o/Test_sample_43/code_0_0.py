from itertools import product

# Define the numbers and target
numbers = [2, 5, 9, 1, 8, 1]
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
    # Construct the expression without parentheses
    expr = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}{ops[4]}{numbers[5]}"
    if safe_eval(expr) == target:
        print(f"Solution without parentheses: {expr}")

    # Try adding parentheses in different positions
    for i in range(1, 6):
        for j in range(i+1, 7):
            expr_with_parens = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}({numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}){ops[4]}{numbers[5]}"
            if safe_eval(expr_with_parens) == target:
                print(f"Solution with parentheses: {expr_with_parens}")