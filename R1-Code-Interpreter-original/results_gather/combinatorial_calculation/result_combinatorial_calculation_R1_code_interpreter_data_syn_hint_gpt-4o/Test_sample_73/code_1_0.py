from itertools import product

# Define the numbers and target
numbers = [7, 3, 4, 8]
target = 24
operations = ['+', '-', '*', '/']

# Function to evaluate an expression safely
def safe_eval(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all possible expressions
def generate_expressions(numbers, operations):
    results = []
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
                results.append(expr)
    return results

# Find all expressions that evaluate to the target
valid_expressions = generate_expressions(numbers, operations)
print(valid_expressions)