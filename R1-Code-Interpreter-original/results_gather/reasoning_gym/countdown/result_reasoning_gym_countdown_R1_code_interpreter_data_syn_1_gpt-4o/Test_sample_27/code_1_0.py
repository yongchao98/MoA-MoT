from itertools import permutations, product

numbers = [39, 12, 33, 86]
target = 778
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    expression = f"({nums[0]} {ops[0]} {nums[1]}) {ops[1]} ({nums[2]} {ops[2]} {nums[3]})"
    try:
        return eval(expression)
    except ZeroDivisionError:
        return None

for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        result = evaluate_expression(num_perm, ops)
        if result == target:
            expression = f"({num_perm[0]} {ops[0]} {num_perm[1]}) {ops[1]} ({num_perm[2]} {ops[2]} {num_perm[3]})"
            print(expression)
            break