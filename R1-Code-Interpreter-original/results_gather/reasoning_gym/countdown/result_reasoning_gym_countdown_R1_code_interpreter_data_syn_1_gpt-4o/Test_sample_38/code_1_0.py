# Remaining numbers
remaining_numbers = [10, 52, 65]

# Check combinations to form 8
from itertools import permutations, product

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression
def evaluate_expression(a, op, b):
    if op == '+':
        return a + b
    elif op == '-':
        return a - b
    elif op == '*':
        return a * b
    elif op == '/':
        return a / b if b != 0 else None

# Try all permutations of numbers and operations
for nums in permutations(remaining_numbers):
    for ops in product(operations, repeat=2):
        # Evaluate (a op1 b) op2 c
        result1 = evaluate_expression(nums[0], ops[0], nums[1])
        if result1 is not None:
            result2 = evaluate_expression(result1, ops[1], nums[2])
            if result2 == 8:
                print(f"Expression: ({nums[0]} {ops[0]} {nums[1]}) {ops[1]} {nums[2]} = 8")