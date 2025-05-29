from itertools import permutations, product

numbers = [82, 38, 56, 49]
target = 588
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
            result /= nums[i]
    return result

for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        if evaluate_expression(num_perm, ops) == target:
            expression = f"{num_perm[0]} {ops[0]} {num_perm[1]} {ops[1]} {num_perm[2]} {ops[2]} {num_perm[3]}"
            print(expression)
            break