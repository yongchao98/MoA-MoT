from itertools import product

# Define the numbers and target
numbers = [7, 4, 5, 7, 5, 7]
target = 100

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression safely
def safe_eval(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Function to generate expressions with parentheses
def generate_expressions(numbers, operations):
    if len(numbers) == 1:
        yield str(numbers[0])
    else:
        for i in range(1, len(numbers)):
            left_numbers = numbers[:i]
            right_numbers = numbers[i:]
            for left_expr in generate_expressions(left_numbers, operations):
                for right_expr in generate_expressions(right_numbers, operations):
                    for op in operations:
                        yield f"({left_expr}{op}{right_expr})"

# Try all combinations of operations
for expr in generate_expressions(numbers, operations):
    if safe_eval(expr) == target:
        print(f"Solution: {expr}")