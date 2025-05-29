from itertools import permutations, product

numbers = [2, 4, 10, 12]
target = 24
operations = ['+', '-', '*', '/']

def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

def generate_expressions(nums, ops):
    # Generate expressions with different groupings
    expressions = [
        f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}",
        f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}",
        f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]})",
        f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]}))",
        f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})"
    ]
    return expressions

for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        for expr in generate_expressions(num_perm, ops):
            if evaluate_expression(expr) == target:
                print(expr)
                break