from itertools import product

# Given numbers
numbers = [9, 2, 8, 9]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression with given operations and parentheses
def evaluate_expression(ops):
    try:
        # Try different parenthesis placements
        expr1 = f"({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} ({numbers[2]} {ops[2]} {numbers[3]})"
        expr2 = f"(({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]}) {ops[2]} {numbers[3]}"
        expr3 = f"{numbers[0]} {ops[0]} (({numbers[1]} {ops[1]} {numbers[2]}) {ops[2]} {numbers[3]})"
        expr4 = f"{numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} ({numbers[2]} {ops[2]} {numbers[3]}))"
        expr5 = f"({numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} {numbers[2]})) {ops[2]} {numbers[3]}"
        
        # Evaluate each expression
        for expr in [expr1, expr2, expr3, expr4, expr5]:
            if eval(expr) == 24:
                return expr
    except ZeroDivisionError:
        pass
    return None

# Try all combinations of operations
for ops in product(operations, repeat=3):
    result = evaluate_expression(ops)
    if result:
        print(result)
        break