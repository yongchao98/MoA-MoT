from itertools import product

numbers = [4, 2, 5]
operations = ['+', '-', '*', '/']
parentheses_positions = [
    (0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2),  # No parentheses
    (0, 3), (1, 3), (2, 3), (3, 3),  # Parentheses around first two
    (0, 4), (1, 4), (2, 4), (3, 4),  # Parentheses around last two
    (0, 5), (1, 5), (2, 5), (3, 5),  # Parentheses around all
]

def evaluate_expression(nums, ops, parens):
    expr = f"{nums[0]}{ops[0]}{nums[1]}{ops[1]}{nums[2]}"
    if parens == (0, 3):
        expr = f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}"
    elif parens == (1, 3):
        expr = f"({nums[0]}{ops[0]}{nums[1]}{ops[1]}{nums[2]})"
    elif parens == (2, 3):
        expr = f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})"
    elif parens == (3, 3):
        expr = f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]})"
    elif parens == (0, 4):
        expr = f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}"
    elif parens == (1, 4):
        expr = f"({nums[0]}{ops[0]}{nums[1]}{ops[1]}{nums[2]})"
    elif parens == (2, 4):
        expr = f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})"
    elif parens == (3, 4):
        expr = f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]})"
    elif parens == (0, 5):
        expr = f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}"
    elif parens == (1, 5):
        expr = f"({nums[0]}{ops[0]}{nums[1]}{ops[1]}{nums[2]})"
    elif parens == (2, 5):
        expr = f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})"
    elif parens == (3, 5):
        expr = f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]})"
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

for ops in product(operations, repeat=2):
    for parens in parentheses_positions:
        if evaluate_expression(numbers, ops, parens) == 10:
            print([numbers[0], ops[0], numbers[1], ops[1], numbers[2]])
            break