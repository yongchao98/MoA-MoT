import itertools

numbers = [62, 37, 46, 53, 30, 41]
target = 157
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

for num_perm in itertools.permutations(numbers):
    for ops in itertools.product(operations, repeat=len(numbers)-1):
        result = evaluate_expression(num_perm, ops)
        if result:
            print(result)
            break