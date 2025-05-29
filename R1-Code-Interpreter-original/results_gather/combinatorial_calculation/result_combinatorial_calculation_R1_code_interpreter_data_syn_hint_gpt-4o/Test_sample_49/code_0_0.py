from itertools import product

# Define the numbers and target
numbers = [8, 9, 3, 1, 8, 5]
target = 100

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression safely
def safe_eval(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all combinations of operations and parentheses
for ops in product(operations, repeat=5):
    # Generate possible expressions with parentheses
    expressions = [
        f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}{numbers[3]}){ops[3]}({numbers[4]}{ops[4]}{numbers[5]})",
        f"(({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}{numbers[2]}){ops[2]}({numbers[3]}{ops[3]}{numbers[4]}{ops[4]}{numbers[5]})",
        f"({numbers[0]}{ops[0]}({numbers[1]}{ops[1]}{numbers[2]})){ops[2]}({numbers[3]}{ops[3]}{numbers[4]}{ops[4]}{numbers[5]})",
        f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}({numbers[3]}{ops[3]}{numbers[4]})){ops[4]}{numbers[5]}",
        f"{numbers[0]}{ops[0]}(({numbers[1]}{ops[1]}{numbers[2]}){ops[2]}({numbers[3]}{ops[3]}{numbers[4]})){ops[4]}{numbers[5]}"
    ]
    
    for expr in expressions:
        if safe_eval(expr) == target:
            print(expr)
            break