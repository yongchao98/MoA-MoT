import itertools

# Define the numbers and target
numbers = [78, 54, 85, 11, 8]
target = 667

# Define the operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(nums, ops):
    try:
        result = nums[0]
        for i in range(4):
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
for num_perm in itertools.permutations(numbers):
    for ops in itertools.product(operations, repeat=4):
        if evaluate_expression(num_perm, ops) == target:
            expression = f"{num_perm[0]} {ops[0]} {num_perm[1]} {ops[1]} {num_perm[2]} {ops[2]} {num_perm[3]} {ops[3]} {num_perm[4]}"
            print(expression)
            break