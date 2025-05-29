from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Evaluate the expression based on the numbers and operations
    result = nums[0]
    expression = str(nums[0])
    for i in range(3):
        if ops[i] == '+':
            result += nums[i + 1]
            expression += f" + {nums[i + 1]}"
        elif ops[i] == '-':
            result -= nums[i + 1]
            expression += f" - {nums[i + 1]}"
        elif ops[i] == '*':
            result *= nums[i + 1]
            expression += f" * {nums[i + 1]}"
        elif ops[i] == '/':
            if nums[i + 1] == 0:
                return None, None  # Avoid division by zero
            result /= nums[i + 1]
            expression += f" / {nums[i + 1]}"
    return result, expression

def find_equation_to_24(numbers):
    operations = ['+', '-', '*', '/']
    for num_perm in permutations(numbers):
        for ops in product(operations, repeat=3):
            result, expression = evaluate_expression(num_perm, ops)
            if result == 24:
                return f"<<<{expression} = 24>>>"
    return "No solution found."

# Input numbers
numbers = [1, 3, 13, 13]

# Find and print the equation
print(find_equation_to_24(numbers))