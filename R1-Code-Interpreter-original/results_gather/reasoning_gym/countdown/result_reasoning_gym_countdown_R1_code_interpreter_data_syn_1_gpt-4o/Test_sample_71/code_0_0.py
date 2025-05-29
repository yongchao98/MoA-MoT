from itertools import permutations, product

numbers = [97, 25, 25, 28, 68, 59]
target = 134
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    result = nums[0]
    for i in range(len(ops)):
        if ops[i] == '+':
            result += nums[i + 1]
        elif ops[i] == '-':
            result -= nums[i + 1]
        elif ops[i] == '*':
            result *= nums[i + 1]
        elif ops[i] == '/':
            if nums[i + 1] == 0:
                return None
            result /= nums[i + 1]
    return result

for num_perm in permutations(numbers):
    for ops in product(operations, repeat=len(numbers)-1):
        if evaluate_expression(num_perm, ops) == target:
            expression = f"{num_perm[0]}"
            for i in range(len(ops)):
                expression += f" {ops[i]} {num_perm[i+1]}"
            print(expression)
            break