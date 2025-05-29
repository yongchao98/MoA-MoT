from itertools import product

# Define the numbers and target
numbers = [9, 7, 6]
target = 10

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all possible expressions with operations and parentheses
def generate_expressions(numbers, operations):
    expressions = []
    for ops in product(operations, repeat=2):
        # Without parentheses
        expr1 = f"{numbers[0]} {ops[0]} {numbers[1]} {ops[1]} {numbers[2]}"
        expressions.append(expr1)
        
        # With parentheses
        expr2 = f"({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]}"
        expr3 = f"{numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} {numbers[2]})"
        expressions.extend([expr2, expr3])
    return expressions

# Find the expression that evaluates to the target
expressions = generate_expressions(numbers, operations)
for expr in expressions:
    if evaluate_expression(expr) == target:
        print(expr)
        break