import itertools

numbers = [65, 95, 66, 55, 27, 15]
target = 137
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    expression = f"{nums[0]}"
    for i in range(1, len(nums)):
        expression += f" {ops[i-1]} {nums[i]}"
    try:
        return eval(expression)
    except ZeroDivisionError:
        return None

def find_expression(numbers, target):
    for nums in itertools.permutations(numbers):
        for ops in itertools.product(operations, repeat=len(numbers)-1):
            if evaluate_expression(nums, ops) == target:
                return f"{nums[0]} {ops[0]} {nums[1]} {ops[1]} {nums[2]} {ops[2]} {nums[3]} {ops[3]} {nums[4]} {ops[4]} {nums[5]}"
    return None

expression = find_expression(numbers, target)
print(expression)