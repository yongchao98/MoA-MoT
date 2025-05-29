from itertools import permutations, product

numbers = [10, 78, 4, 86, 92]
target = 340
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    expression = f"{nums[0]}"
    for i in range(4):
        expression += f" {ops[i]} {nums[i+1]}"
    try:
        return eval(expression)
    except ZeroDivisionError:
        return None

for num_perm in permutations(numbers):
    for ops in product(operations, repeat=4):
        if evaluate_expression(num_perm, ops) == target:
            expression = f"{num_perm[0]}"
            for i in range(4):
                expression += f" {ops[i]} {num_perm[i+1]}"
            print(expression)
            break