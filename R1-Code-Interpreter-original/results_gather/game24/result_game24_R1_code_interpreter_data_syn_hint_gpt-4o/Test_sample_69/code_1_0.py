import itertools
from sympy import symbols, Eq, solve

# Define the numbers and target
numbers = [4, 5, 7, 7]
target = 24

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(nums, ops):
    x, y, z, w = symbols('x y z w')
    expr = f"(({x} {ops[0]} {y}) {ops[1]} ({z} {ops[2]} {w}))"
    expr = expr.replace('x', str(nums[0])).replace('y', str(nums[1])).replace('z', str(nums[2])).replace('w', str(nums[3]))
    try:
        result = eval(expr)
        return result, expr
    except ZeroDivisionError:
        return None, None

# Try all permutations of numbers and operations
for num_perm in itertools.permutations(numbers):
    for ops in itertools.product(operations, repeat=3):
        result, expression = evaluate_expression(num_perm, ops)
        if result == target:
            print(expression)
            break