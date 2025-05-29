from itertools import permutations

numbers = [17, 43, 31, 5, 94]
target = 651

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
            if nums[i + 1] == 0:
                return None
            result /= nums[i + 1]
    return result

# Try all permutations of numbers and operations
for num_perm in permutations(numbers):
    for ops in permutations(['+', '-', '*', '/'], 4):
        if evaluate_expression(num_perm, ops) == target:
            expression = f"{num_perm[0]}"
            for i in range(4):
                expression += f" {ops[i]} {num_perm[i + 1]}"
            print(expression)
            break