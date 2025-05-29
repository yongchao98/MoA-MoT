from sympy import symbols, Eq, solve
from itertools import product

# Define the numbers and operations
numbers = [6, 9, 6, 7, 2, 8]
operations = ['+', '-', '*', '/']

# Define a function to evaluate expressions using SymPy
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all combinations of operations
for ops in product(operations, repeat=5):
    # Generate possible expressions with parentheses
    expressions = [
        f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}{numbers[3]}){ops[3]}({numbers[4]}{ops[4]}{numbers[5]})",
        f"(({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}{numbers[2]}){ops[2]}({numbers[3]}{ops[3]}{numbers[4]}){ops[4]}{numbers[5]}",
        f"({numbers[0]}{ops[0]}({numbers[1]}{ops[1]}{numbers[2]})){ops[2]}({numbers[3]}{ops[3]}{numbers[4]}){ops[4]}{numbers[5]}",
        f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}({numbers[3]}{ops[3]}{numbers[4]})){ops[4]}{numbers[5]}",
        f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}{numbers[3]}){ops[3]}({numbers[4]}{ops[4]}{numbers[5]})",
        f"{numbers[0]}{ops[0]}(({numbers[1]}{ops[1]}{numbers[2]}){ops[2]}{numbers[3]}){ops[3]}({numbers[4]}{ops[4]}{numbers[5]})",
        f"{numbers[0]}{ops[0]}({numbers[1]}{ops[1]}({numbers[2]}{ops[2]}{numbers[3]})){ops[3]}({numbers[4]}{ops[4]}{numbers[5]})",
        f"{numbers[0]}{ops[0]}({numbers[1]}{ops[1]}{numbers[2]}){ops[2]}(({numbers[3]}{ops[3]}{numbers[4]}){ops[4]}{numbers[5]})",
    ]
    
    for expr in expressions:
        if evaluate_expression(expr) == 100:
            print(expr)
            break