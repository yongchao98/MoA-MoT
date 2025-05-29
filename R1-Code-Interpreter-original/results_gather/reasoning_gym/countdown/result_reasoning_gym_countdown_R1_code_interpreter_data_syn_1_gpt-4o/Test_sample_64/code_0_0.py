from itertools import permutations, product

numbers = [57, 8, 16, 22, 65]
target = 216
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    result = nums[0]
    for i in range(1, len(nums)):
        if ops[i-1] == '+':
            result += nums[i]
        elif ops[i-1] == '-':
            result -= nums[i]
        elif ops[i-1] == '*':
            result *= nums[i]
        elif ops[i-1] == '/':
            if nums[i] == 0 or result % nums[i] != 0:
                return None
            result //= nums[i]
    return result

for num_perm in permutations(numbers):
    for ops in product(operations, repeat=4):
        if evaluate_expression(num_perm, ops) == target:
            expression = f"{num_perm[0]}"
            for i in range(4):
                expression += f" {ops[i]} {num_perm[i+1]}"
            print(expression)
            break