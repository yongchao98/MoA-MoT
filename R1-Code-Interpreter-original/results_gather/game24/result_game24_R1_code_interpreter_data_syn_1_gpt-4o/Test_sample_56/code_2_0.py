from itertools import permutations, product

numbers = [1, 9, 12, 12]
target = 24
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    try:
        # Evaluate the expression with the given numbers and operations
        result = nums[0]
        for i in range(3):
            if ops[i] == '+':
                result += nums[i+1]
            elif ops[i] == '-':
                result -= nums[i+1]
            elif ops[i] == '*':
                result *= nums[i+1]
            elif ops[i] == '/':
                result /= nums[i+1]
        return result
    except ZeroDivisionError:
        return None

def format_expression(nums, ops):
    # Format the expression with parentheses to ensure correct order of operations
    return f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})"

for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        if evaluate_expression(num_perm, ops) == target:
            expression = format_expression(num_perm, ops)
            print(expression)
            break