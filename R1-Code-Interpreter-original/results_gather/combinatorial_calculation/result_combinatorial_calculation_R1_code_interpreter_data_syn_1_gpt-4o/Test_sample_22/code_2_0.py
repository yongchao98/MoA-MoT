from itertools import product

# Given numbers
numbers = [8, 9, 7, 2, 9]

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate expression with parentheses
def evaluate_expression(nums, ops, parens):
    expr = f"{parens[0]}{nums[0]}{ops[0]}{parens[1]}{nums[1]}{parens[2]}{ops[1]}{parens[3]}{nums[2]}{parens[4]}{ops[2]}{parens[5]}{nums[3]}{parens[6]}{ops[3]}{parens[7]}{nums[4]}{parens[8]}"
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None
    except SyntaxError:
        return None

# Generate all combinations of operations and parentheses
for ops in product(operations, repeat=4):
    for parens in product(['', '(', ')'], repeat=9):
        if parens.count('(') == parens.count(')'):
            result = evaluate_expression(numbers, ops, parens)
            if result == 78:
                print([parens[0], numbers[0], ops[0], parens[1], numbers[1], parens[2], ops[1], parens[3], numbers[2], parens[4], ops[2], parens[5], numbers[3], parens[6], ops[3], parens[7], numbers[4], parens[8]])
                break