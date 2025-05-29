from itertools import permutations, product

numbers = [89, 75, 84, 56, 94, 73]
target = 133
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
            if nums[i + 1] == 0 or result % nums[i + 1] != 0:
                return None
            result //= nums[i + 1]
    return result

def find_expression(numbers, target):
    for num_perm in permutations(numbers):
        for ops in product(operations, repeat=len(numbers) - 1):
            if evaluate_expression(num_perm, ops) == target:
                expression = f"{num_perm[0]}"
                for i in range(len(ops)):
                    expression += f" {ops[i]} {num_perm[i + 1]}"
                return expression
    return None

expression = find_expression(numbers, target)
print(expression)