# Let's try different combinations of numbers and operations to reach 537
from itertools import permutations, product

numbers = [51, 46, 36, 38, 21]
target = 537

# Function to evaluate an expression
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
            result /= nums[i + 1]
    return result

# Try all permutations of numbers and operations
for num_perm in permutations(numbers):
    for ops in product('+-*/', repeat=4):
        try:
            if evaluate_expression(num_perm, ops) == target:
                expression = f"{num_perm[0]} {ops[0]} {num_perm[1]} {ops[1]} {num_perm[2]} {ops[2]} {num_perm[3]} {ops[3]} {num_perm[4]}"
                print(expression)
        except ZeroDivisionError:
            continue