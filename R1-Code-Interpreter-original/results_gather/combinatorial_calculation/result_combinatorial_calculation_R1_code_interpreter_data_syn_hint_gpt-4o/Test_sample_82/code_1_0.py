from itertools import permutations, product

# Define the numbers and target
numbers = [3, 1, 9, 6]
target = 24
operations = ['+', '-', '*', '/']

# Function to evaluate expression safely
def safe_eval(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all possible expressions with parentheses
def generate_expressions(numbers, operations):
    expressions = []
    for ops in product(operations, repeat=3):
        # Without parentheses
        expressions.append(f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}")
        # With different parentheses
        expressions.append(f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}{numbers[2]}{ops[2]}{numbers[3]}")
        expressions.append(f"{numbers[0]}{ops[0]}({numbers[1]}{ops[1]}{numbers[2]}){ops[2]}{numbers[3]}")
        expressions.append(f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}({numbers[2]}{ops[2]}{numbers[3]})")
        expressions.append(f"({numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}){ops[2]}{numbers[3]}")
        expressions.append(f"{numbers[0]}{ops[0]}({numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]})")
        expressions.append(f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}{numbers[3]})")
    return expressions

# Check all expressions
for expr in generate_expressions(numbers, operations):
    if safe_eval(expr) == target:
        print(expr)
        break