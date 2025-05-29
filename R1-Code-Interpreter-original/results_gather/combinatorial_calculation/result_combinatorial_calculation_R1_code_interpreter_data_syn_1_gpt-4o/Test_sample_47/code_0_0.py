from itertools import product

# Define the numbers and target
numbers = [7, 4, 5, 5, 5, 5]
target = 100

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression safely
def safe_eval(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all combinations of operations
for ops in product(operations, repeat=5):
    # Generate expression without parentheses
    expr = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}{ops[4]}{numbers[5]}"
    if safe_eval(expr) == target:
        print(f"Expression: {expr}")
        break

    # Generate expressions with parentheses
    for i in range(5):
        for j in range(i+2, 7):
            expr = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}{ops[4]}{numbers[5]}"
            expr = expr[:i] + '(' + expr[i:j] + ')' + expr[j:]
            if safe_eval(expr) == target:
                print(f"Expression: {expr}")
                break