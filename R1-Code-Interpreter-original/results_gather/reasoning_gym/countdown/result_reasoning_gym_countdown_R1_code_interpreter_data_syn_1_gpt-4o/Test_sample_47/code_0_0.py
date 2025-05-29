from itertools import permutations, product

# Define the numbers and target
numbers = [5, 23, 35, 9]
target = 809
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(nums, ops):
    try:
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
    except ZeroDivisionError:
        return None

# Try all permutations of numbers and operations
for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        if evaluate_expression(num_perm, ops) == target:
            expression = f"{num_perm[0]} {ops[0]} {num_perm[1]} {ops[1]} {num_perm[2]} {ops[2]} {num_perm[3]}"
            print(expression)
            break