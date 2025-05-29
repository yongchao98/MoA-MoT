from itertools import permutations, product

numbers = [75, 23, 29, 32, 93, 75]
target = 281
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
        if evaluate_expression(num_perm, ops) == target:
            expression = str(num_perm[0])
            for i in range(1, len(num_perm)):
                expression += ops[i-1] + str(num_perm[i])
            print(expression)
            break