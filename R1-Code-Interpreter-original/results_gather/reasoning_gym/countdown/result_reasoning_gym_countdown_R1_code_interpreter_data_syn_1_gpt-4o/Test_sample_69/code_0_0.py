from itertools import permutations, product

numbers = [61, 94, 88, 84, 47, 73]
target = 165
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    expression = str(nums[0])
    for i in range(len(ops)):
        expression += ops[i] + str(nums[i+1])
    try:
        return eval(expression)
    except ZeroDivisionError:
        return None

for num_perm in permutations(numbers):
    for ops in product(operations, repeat=len(numbers)-1):
        if evaluate_expression(num_perm, ops) == target:
            expression = str(num_perm[0])
            for i in range(len(ops)):
                expression += ops[i] + str(num_perm[i+1])
            print(expression)
            break