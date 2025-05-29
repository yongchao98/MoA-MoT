from itertools import permutations, product

numbers = [1, 6, 6, 11]
target = 24
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    try:
        # Try different parenthesized expressions
        expressions = [
            f"(({nums[0]} {ops[0]} {nums[1]}) {ops[1]} ({nums[2]} {ops[2]} {nums[3]}))",
            f"(({nums[0]} {ops[0]} ({nums[1]} {ops[1]} {nums[2]})) {ops[2]} {nums[3]})",
            f"({nums[0]} {ops[0]} (({nums[1]} {ops[1]} {nums[2]}) {ops[2]} {nums[3]}))",
            f"({nums[0]} {ops[0]} ({nums[1]} {ops[1]} ({nums[2]} {ops[2]} {nums[3]})))"
        ]
        for expr in expressions:
            if eval(expr) == target:
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