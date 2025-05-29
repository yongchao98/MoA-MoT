from itertools import permutations, product

numbers = [2, 6, 7, 9]
target = 24
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    try:
        # Construct the expression with parentheses
        expr1 = f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]}))"
        expr2 = f"(({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]})"
        expr3 = f"({nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]}))"
        expr4 = f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]})))"
        expr5 = f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}"
        
        for expr in [expr1, expr2, expr3, expr4, expr5]:
            if abs(eval(expr) - target) < 1e-6:
                return expr
    except ZeroDivisionError:
        return None
    return None

for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        result = evaluate_expression(num_perm, ops)
        if result:
            print(result)
            break