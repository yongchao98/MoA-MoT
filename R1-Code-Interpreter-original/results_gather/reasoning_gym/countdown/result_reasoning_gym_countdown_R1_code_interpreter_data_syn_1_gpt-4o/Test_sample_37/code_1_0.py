from itertools import permutations, product

numbers = [22, 78, 13, 45, 5, 9]
target = 264
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    expression = f"{nums[0]}"
    for i in range(len(ops)):
        expression += f" {ops[i]} {nums[i+1]}"
    try:
        return eval(expression)
    except ZeroDivisionError:
        return None

for num_perm in permutations(numbers):
    for ops in product(operations, repeat=5):
        if evaluate_expression(num_perm, ops) == target:
            expression = f"{num_perm[0]}"
            for i in range(len(ops)):
                expression += f" {ops[i]} {num_perm[i+1]}"
            print(expression)
            break