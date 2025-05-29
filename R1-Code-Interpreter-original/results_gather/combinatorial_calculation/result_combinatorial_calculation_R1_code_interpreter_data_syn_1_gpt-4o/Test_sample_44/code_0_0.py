from itertools import product

# Define the numbers and target
numbers = [9, 9, 4, 6, 2, 8]
target = 100

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression with error handling for division by zero
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all combinations of operations
for ops in product(operations, repeat=5):
    # Generate possible expressions with parentheses
    expressions = [
        f"(({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}{numbers[2]}){ops[2]}({numbers[3]}{ops[3]}{numbers[4]}){ops[4]}{numbers[5]}",
        f"({numbers[0]}{ops[0]}({numbers[1]}{ops[1]}{numbers[2]})){ops[2]}({numbers[3]}{ops[3]}{numbers[4]}){ops[4]}{numbers[5]}",
        f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}({numbers[3]}{ops[3]}{numbers[4]})){ops[4]}{numbers[5]}",
        f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}{numbers[3]}){ops[3]}({numbers[4]}{ops[4]}{numbers[5]})",
        f"(({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}{numbers[2]}){ops[2]}{numbers[3]}{ops[3]}({numbers[4]}{ops[4]}{numbers[5]})",
        f"{numbers[0]}{ops[0]}(({numbers[1]}{ops[1]}{numbers[2]}){ops[2]}{numbers[3]}){ops[3]}({numbers[4]}{ops[4]}{numbers[5]})",
        f"{numbers[0]}{ops[0]}({numbers[1]}{ops[1]}({numbers[2]}{ops[2]}{numbers[3]})){ops[3]}({numbers[4]}{ops[4]}{numbers[5]})",
        f"{numbers[0]}{ops[0]}({numbers[1]}{ops[1]}{numbers[2]}){ops[2]}(({numbers[3]}{ops[3]}{numbers[4]}){ops[4]}{numbers[5]})",
    ]
    
    # Evaluate each expression
    for expr in expressions:
        if evaluate_expression(expr) == target:
            print(expr)
            break