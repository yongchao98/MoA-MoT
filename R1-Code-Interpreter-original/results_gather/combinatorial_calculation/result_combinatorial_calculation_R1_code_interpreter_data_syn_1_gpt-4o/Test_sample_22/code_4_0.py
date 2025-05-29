from itertools import product

# Given numbers
numbers = [8, 9, 7, 2, 9]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression with parentheses
def evaluate_expression(nums, ops, parens):
    expr = f"{parens[0]}{nums[0]}{parens[1]}{ops[0]}{parens[2]}{nums[1]}{parens[3]}{ops[1]}{parens[4]}{nums[2]}{parens[5]}{ops[2]}{parens[6]}{nums[3]}{parens[7]}{ops[3]}{parens[8]}{nums[4]}{parens[9]}"
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None
    except SyntaxError:
        return None

# Generate all combinations of operations and parentheses
for ops in product(operations, repeat=4):
    for parens in product(['', '(', ')'], repeat=10):
        if parens.count('(') == parens.count(')'):
            result = evaluate_expression(numbers, ops, parens)
            if result == 78:
                print([parens[0], numbers[0], parens[1], ops[0], parens[2], numbers[1], parens[3], ops[1], parens[4], numbers[2], parens[5], ops[2], parens[6], numbers[3], parens[7], ops[3], parens[8], numbers[4], parens[9]])
                break