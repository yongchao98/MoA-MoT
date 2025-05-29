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

# Generate all possible combinations of operations
for ops in product(operations, repeat=5):
    # Try different parenthesis placements
    expressions = [
        f"(({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]}) {ops[2]} ({numbers[3]} {ops[3]} {numbers[4]}) {ops[4]} {numbers[5]}",
        f"({numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} {numbers[2]})) {ops[2]} ({numbers[3]} {ops[3]} {numbers[4]}) {ops[4]} {numbers[5]}",
        f"({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} ({numbers[2]} {ops[2]} ({numbers[3]} {ops[3]} {numbers[4]})) {ops[4]} {numbers[5]}",
        f"({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]} {ops[2]} ({numbers[3]} {ops[3]} ({numbers[4]} {ops[4]} {numbers[5]}))",
        f"(({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]}) {ops[2]} {numbers[3]} {ops[3]} ({numbers[4]} {ops[4]} {numbers[5]})",
    ]
    
    for expr in expressions:
        if safe_eval(expr) == target:
            print(expr)
            break