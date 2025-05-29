from sympy import symbols, Eq, solve
from itertools import product

# Define the numbers as symbols
a, b, c, d, e, f = symbols('a b c d e f')
numbers = [9, 6, 1, 7, 1, 3]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to generate expressions
def generate_expressions(nums, ops):
    exprs = [
        f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]}){ops[3]}({nums[4]}{ops[4]}{nums[5]})",
        f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}({nums[3]}{ops[3]}{nums[4]}{ops[4]}{nums[5]})",
        f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}({nums[3]}{ops[3]}{nums[4]}{ops[4]}{nums[5]})",
        f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}({nums[3]}{ops[3]}{nums[4]})){ops[4]}{nums[5]}",
        f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]}){ops[3]}({nums[4]}{ops[4]}{nums[5]})",
        f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]})){ops[3]}({nums[4]}{ops[4]}{nums[5]})"
    ]
    return exprs

# Generate all combinations of operations
for ops in product(operations, repeat=5):
    expressions = generate_expressions(numbers, ops)
    
    for expr in expressions:
        try:
            if eval(expr) == 100:
                print(expr)
                break
        except ZeroDivisionError:
            continue