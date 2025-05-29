from itertools import permutations, product

numbers = [70, 74, 67, 86, 80, 19]
target = 358
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    expression = str(nums[0])
    for i in range(1, len(nums)):
        expression += ops[i-1] + str(nums[i])
    try:
        if eval(expression) == target:
            return expression
    except ZeroDivisionError:
        pass
    return None

for num_perm in permutations(numbers):
    for ops in product(operations, repeat=len(numbers)-1):
        result = evaluate_expression(num_perm, ops)
        if result:
            print(result)
            break