from itertools import product

# Given numbers
numbers = [8, 9, 7, 2, 9]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression with parentheses
def evaluate_expression(nums, ops, parens):
    expr = f"{parens[0]}{nums[0]}{ops[0]}{nums[1]}{parens[1]}{ops[1]}{parens[2]}{nums[2]}{parens[3]}{ops[2]}{parens[4]}{nums[3]}{parens[5]}{ops[3]}{nums[4]}{parens[6]}"
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all combinations of operations and parentheses
for ops in product(operations, repeat=4):
    for parens in product(['', '(', ')'], repeat=7):
        if parens.count('(') == parens.count(')'):
            result = evaluate_expression(numbers, ops, parens)
            if result == 78:
                print([numbers[0], ops[0], parens[0], numbers[1], parens[1], ops[1], parens[2], numbers[2], parens[3], ops[2], parens[4], numbers[3], parens[5], ops[3], numbers[4], parens[6]])
                break