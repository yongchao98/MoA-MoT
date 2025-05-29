from itertools import permutations
from operator import add, sub, mul, truediv

def evaluate_expression(nums, ops):
    # Evaluate the expression based on the current order of numbers and operations
    result = nums[0]
    for i in range(3):
        if ops[i] == '+':
            result = add(result, nums[i + 1])
        elif ops[i] == '-':
            result = sub(result, nums[i + 1])
        elif ops[i] == '*':
            result = mul(result, nums[i + 1])
        elif ops[i] == '/':
            if nums[i + 1] == 0:
                return None
            result = truediv(result, nums[i + 1])
    return result

def find_expression(nums):
    operations = ['+', '-', '*', '/']
    for num_perm in permutations(nums):
        for ops in permutations(operations, 3):
            if evaluate_expression(num_perm, ops) == 24:
                return f"<<<({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}{num_perm[2]}{ops[2]}{num_perm[3]}>>>"
    return "No solution found"

nums = [2, 4, 10, 12]
print(find_expression(nums))