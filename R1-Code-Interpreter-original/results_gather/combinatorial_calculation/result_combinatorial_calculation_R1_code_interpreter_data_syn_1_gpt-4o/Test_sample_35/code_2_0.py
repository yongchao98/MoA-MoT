from itertools import product

# Given numbers
numbers = [8, 5, 3, 2, 5]

# Possible operations
operations = ['+', '-', '*', '/']

# Target value
target = 69

# Function to evaluate expression with given operations and parentheses
def evaluate_expression(numbers, ops, parens):
    expr = f"{parens[0]}{numbers[0]}{ops[0]}{parens[1]}{numbers[1]}{parens[2]}{ops[1]}{numbers[2]}{parens[3]}{ops[2]}{numbers[3]}{parens[4]}{ops[3]}{numbers[4]}{parens[5]}"
    try:
        return eval(expr)
    except (ZeroDivisionError, SyntaxError, TypeError):
        return None

# Try all combinations of operations and parentheses
for ops in product(operations, repeat=4):
    for parens in product(['', '(', ')'], repeat=6):
        # Ensure parentheses are balanced and correctly placed
        if parens.count('(') == parens.count(')') and not any(parens[i] == ')' and parens[i-1] == '(' for i in range(1, 6)):
            result = evaluate_expression(numbers, ops, parens)
            if result == target:
                print([parens[0], numbers[0], ops[0], parens[1], numbers[1], parens[2], ops[1], numbers[2], parens[3], ops[2], numbers[3], parens[4], ops[3], numbers[4], parens[5]])
                break