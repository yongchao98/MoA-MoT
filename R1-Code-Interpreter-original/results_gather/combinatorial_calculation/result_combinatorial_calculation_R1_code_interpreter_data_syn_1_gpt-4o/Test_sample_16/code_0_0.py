from itertools import product

# Define the numbers and target
numbers = [5, 1, 3, 2]
target = 24

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(numbers, ops):
    try:
        # Construct the expression with the given operations
        expression = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}"
        # Evaluate the expression
        return eval(expression)
    except ZeroDivisionError:
        return None

# Try all combinations of operations
for ops in product(operations, repeat=3):
    # Evaluate without parentheses
    if evaluate_expression(numbers, ops) == target:
        print([numbers[0], ops[0], numbers[1], ops[1], numbers[2], ops[2], numbers[3]])
    
    # Evaluate with different parentheses placements
    # (a op1 b) op2 (c op3 d)
    if eval(f"({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}({numbers[2]}{ops[2]}{numbers[3]})") == target:
        print(['(', numbers[0], ops[0], numbers[1], ')', ops[1], '(', numbers[2], ops[2], numbers[3], ')'])
    
    # ((a op1 b) op2 c) op3 d
    if eval(f"(({numbers[0]}{ops[0]}{numbers[1]}){ops[1]}{numbers[2]}){ops[2]}{numbers[3]}") == target:
        print(['(', '(', numbers[0], ops[0], numbers[1], ')', ops[1], numbers[2], ')', ops[2], numbers[3]])
    
    # (a op1 (b op2 c)) op3 d
    if eval(f"({numbers[0]}{ops[0]}({numbers[1]}{ops[1]}{numbers[2]})){ops[2]}{numbers[3]}") == target:
        print(['(', numbers[0], ops[0], '(', numbers[1], ops[1], numbers[2], ')', ')', ops[2], numbers[3]])
    
    # a op1 ((b op2 c) op3 d)
    if eval(f"{numbers[0]}{ops[0]}(({numbers[1]}{ops[1]}{numbers[2]}){ops[2]}{numbers[3]})") == target:
        print([numbers[0], ops[0], '(', '(', numbers[1], ops[1], numbers[2], ')', ops[2], numbers[3], ')'])
