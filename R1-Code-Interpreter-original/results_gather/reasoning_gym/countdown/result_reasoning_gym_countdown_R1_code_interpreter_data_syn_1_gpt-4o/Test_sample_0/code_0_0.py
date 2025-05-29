from itertools import permutations, product

# Define the numbers and target
numbers = [11, 26, 43, 58]
target = 138

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(nums, ops):
    result = nums[0]
    for i in range(3):
        if ops[i] == '+':
            result += nums[i+1]
        elif ops[i] == '-':
            result -= nums[i+1]
        elif ops[i] == '*':
            result *= nums[i+1]
        elif ops[i] == '/':
            if nums[i+1] == 0:
                return None
            result /= nums[i+1]
    return result

# Try all permutations of numbers and operations
for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        if evaluate_expression(num_perm, ops) == target:
            expression = f"{num_perm[0]} {ops[0]} {num_perm[1]} {ops[1]} {num_perm[2]} {ops[2]} {num_perm[3]}"
            print(expression)
            break