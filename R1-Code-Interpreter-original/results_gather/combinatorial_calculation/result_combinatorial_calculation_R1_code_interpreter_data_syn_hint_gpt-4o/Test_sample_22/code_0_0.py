from itertools import product

# Define the numbers and target
numbers = [8, 9, 7, 2, 9]
target = 78
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all possible combinations of operations
for ops in product(operations, repeat=4):
    # Generate possible expressions with parentheses
    expressions = [
        f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}{numbers[3]}){ops[3]}{numbers[4]}",
        f"(({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}{numbers[2]}){ops[2]}({numbers[3]}{ops[3]}{numbers[4]})",
        f"({numbers[0]}{ops[0]}({numbers[1]}{ops[1]}{numbers[2]})){ops[2]}({numbers[3]}{ops[3]}{numbers[4]})",
        f"{numbers[0]}{ops[0]}(({numbers[1]}{ops[1]}{numbers[2]}){ops[2]}{numbers[3]}){ops[3]}{numbers[4]}",
        f"{numbers[0]}{ops[0]}({numbers[1]}{ops[1]}({numbers[2]}{ops[2]}{numbers[3]})){ops[3]}{numbers[4]}"
    ]
    
    # Evaluate each expression
    for expr in expressions:
        if evaluate_expression(expr) == target:
            print(expr)
            break