from itertools import product

# Given numbers
numbers = [8, 3, 7, 1, 9]
# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression safely
def safe_eval(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all combinations of operations
for ops in product(operations, repeat=4):
    # Generate expressions with different parenthesis placements
    expressions = [
        f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}{numbers[3]}){ops[3]}{numbers[4]}",
        f"(({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}{numbers[2]}){ops[2]}({numbers[3]}{ops[3]}{numbers[4]})",
        f"({numbers[0]}{ops[0]}({numbers[1]}{ops[1]}{numbers[2]})){ops[2]}({numbers[3]}{ops[3]}{numbers[4]})",
        f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}({numbers[3]}{ops[3]}{numbers[4]}))",
        f"{numbers[0]}{ops[0]}(({numbers[1]}{ops[1]}{numbers[2]}){ops[2]}({numbers[3]}{ops[3]}{numbers[4]}))",
        f"{numbers[0]}{ops[0]}({numbers[1]}{ops[1]}({numbers[2]}{ops[2]}({numbers[3]}{ops[3]}{numbers[4]})))"
    ]
    
    # Evaluate each expression
    for expr in expressions:
        if safe_eval(expr) == 78:
            print(expr)
            break