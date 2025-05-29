from itertools import permutations, product

numbers = [100, 59, 82, 76, 18, 58]
target = 986
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    expression = f"({nums[0]} {ops[0]} {nums[1]})"
    for i in range(1, len(ops)):
        expression = f"({expression} {ops[i]} {nums[i+1]})"
    try:
        return eval(expression)
    except ZeroDivisionError:
        return None

for num_perm in permutations(numbers):
    for ops in product(operations, repeat=5):
        if evaluate_expression(num_perm, ops) == target:
            expression = f"({num_perm[0]} {ops[0]} {num_perm[1]})"
            for i in range(1, len(ops)):
                expression = f"({expression} {ops[i]} {num_perm[i+1]})"
            print(expression)
            break