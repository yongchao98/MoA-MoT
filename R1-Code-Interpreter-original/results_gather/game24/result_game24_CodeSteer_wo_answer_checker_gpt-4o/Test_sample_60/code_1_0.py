from itertools import permutations, product
from sympy import symbols, Eq, solve

# Define the numbers and target
numbers = [2, 8, 8, 13]
target = 24

# Define the operations
operations = ['+', '-', '*', '/']

# Define a function to evaluate expressions
def evaluate_expression(nums, ops):
    x, y, z, w = symbols('x y z w')
    expr = f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]}))"
    try:
        # Evaluate the expression
        if solve(Eq(eval(expr), target)):
            return expr
    except ZeroDivisionError:
        return None
    return None

# Try all permutations of numbers and combinations of operations
for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        result = evaluate_expression(num_perm, ops)
        if result:
            print(f"<<<{result}>>>")
            break