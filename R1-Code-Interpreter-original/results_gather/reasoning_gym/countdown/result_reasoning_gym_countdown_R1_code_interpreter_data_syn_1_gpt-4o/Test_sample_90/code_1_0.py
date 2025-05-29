from itertools import permutations, product

numbers = [10, 43, 99, 100, 36, 10]
target = 190
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    expression = f"{nums[0]}"
    for i in range(1, len(nums)):
        expression += f" {ops[i-1]} {nums[i]}"
    try:
        return eval(expression)
    except ZeroDivisionError:
        return None

for num_perm in permutations(numbers):
    for ops in product(operations, repeat=5):
        if evaluate_expression(num_perm, ops) == target:
            expression = f"{num_perm[0]}"
            for i in range(1, len(num_perm)):
                expression += f" {ops[i-1]} {num_perm[i]}"
            print(expression)
            break