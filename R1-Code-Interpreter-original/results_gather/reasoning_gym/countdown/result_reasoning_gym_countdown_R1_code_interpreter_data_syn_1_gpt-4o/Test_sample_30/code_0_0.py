from itertools import permutations, product

numbers = [86, 40, 60, 23, 9, 31]
target = 162
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    expression = f"{nums[0]}"
    for i in range(len(ops)):
        expression += f" {ops[i]} {nums[i+1]}"
    try:
        if eval(expression) == target:
            return expression
    except ZeroDivisionError:
        pass
    return None

for nums in permutations(numbers):
    for ops in product(operations, repeat=len(numbers)-1):
        result = evaluate_expression(nums, ops)
        if result:
            print(result)
            break