from itertools import permutations, product

numbers = [92, 32, 88, 33, 8, 10]
target = 935
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    expression = str(nums[0])
    for i in range(len(ops)):
        expression += ops[i] + str(nums[i+1])
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