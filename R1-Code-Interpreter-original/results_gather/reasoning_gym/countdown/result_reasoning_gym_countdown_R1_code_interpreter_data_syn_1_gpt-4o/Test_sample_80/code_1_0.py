from itertools import permutations, product

numbers = [75, 90, 80, 39, 86]
target = 751
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    expression = str(nums[0])
    for i in range(1, len(nums)):
        expression += ops[i-1] + str(nums[i])
    try:
        return eval(expression)
    except ZeroDivisionError:
        return None

for num_perm in permutations(numbers):
    for ops in product(operations, repeat=len(numbers)-1):
        result = evaluate_expression(num_perm, ops)
        if result == target:
            expression = f"{num_perm[0]}"
            for i in range(1, len(num_perm)):
                expression += f" {ops[i-1]} {num_perm[i]}"
            print(expression)
            break