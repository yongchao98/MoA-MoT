from itertools import product

# Define the numbers and target
numbers = [3, 8, 7, 7]
target = 24

# Define the operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression safely
def safe_eval(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Try all combinations of operations and parentheses
def find_expression(numbers, target):
    for ops in product(operations, repeat=3):
        # Generate expressions with different parenthesis placements
        expressions = [
            f"({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} ({numbers[2]} {ops[2]} {numbers[3]})",
            f"(({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]}) {ops[2]} {numbers[3]}",
            f"({numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} {numbers[2]})) {ops[2]} {numbers[3]}",
            f"{numbers[0]} {ops[0]} (({numbers[1]} {ops[1]} {numbers[2]}) {ops[2]} {numbers[3]})",
            f"{numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} ({numbers[2]} {ops[2]} {numbers[3]}))"
        ]
        
        for expr in expressions:
            if safe_eval(expr) == target:
                return expr

# Find and print the expression that evaluates to the target
result = find_expression(numbers, target)
print(result)