from itertools import permutations, product

numbers = [68, 84, 80, 35, 3, 34]
target = 765
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    expression = f"{nums[0]}"
    for i in range(1, len(nums)):
        expression += f" {ops[i-1]} {nums[i]}"
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