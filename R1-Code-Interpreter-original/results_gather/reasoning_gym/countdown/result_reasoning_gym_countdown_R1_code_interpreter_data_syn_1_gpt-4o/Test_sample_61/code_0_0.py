import itertools

numbers = [16, 5, 75, 80, 21]
target = 635
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    expression = f"({nums[0]} {ops[0]} {nums[1]}) {ops[1]} ({nums[2]} {ops[2]} {nums[3]}) {ops[3]} {nums[4]}"
    try:
        if eval(expression) == target:
            return expression
    except ZeroDivisionError:
        pass
    return None

for num_perm in itertools.permutations(numbers):
    for ops in itertools.product(operations, repeat=4):
        result = evaluate_expression(num_perm, ops)
        if result:
            print(result)
            break