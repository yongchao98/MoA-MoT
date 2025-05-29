from itertools import product

# Define the numbers and target
numbers = [7, 7, 9]
target = 10

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression with given operations and parentheses
def evaluate_expression(nums, ops, parens):
    expr = f"{parens[0]}{nums[0]}{ops[0]}{parens[1]}{nums[1]}{parens[2]}{ops[1]}{nums[2]}{parens[3]}"
    try:
        return eval(expr)
    except (ZeroDivisionError, SyntaxError, TypeError):
        return None

# Check all combinations of operations and parentheses
for ops in product(operations, repeat=2):
    for parens in product(['', '(', ')'], repeat=4):
        # Ensure parentheses are valid
        if parens.count('(') == parens.count(')') and (parens.count('(') == 0 or parens.index('(') < parens.index(')')):
            result = evaluate_expression(numbers, ops, parens)
            if result == target:
                print([parens[0], numbers[0], ops[0], parens[1], numbers[1], parens[2], ops[1], numbers[2], parens[3]])
