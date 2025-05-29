from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Evaluate the expression based on the current permutation of numbers and operations
    result = nums[0]
    expression = str(nums[0])
    for i in range(3):
        if ops[i] == '+':
            result += nums[i + 1]
            expression += f"+{nums[i + 1]}"
        elif ops[i] == '-':
            result -= nums[i + 1]
            expression += f"-{nums[i + 1]}"
        elif ops[i] == '*':
            result *= nums[i + 1]
            expression += f"*{nums[i + 1]}"
        elif ops[i] == '/':
            if nums[i + 1] == 0:
                return None, None  # Avoid division by zero
            result /= nums[i + 1]
            expression += f"/{nums[i + 1]}"
    return result, expression

def find_expression_to_24(numbers):
    # Generate all permutations of the numbers
    for nums in permutations(numbers):
        # Generate all combinations of operations
        for ops in product('+-*/', repeat=3):
            result, expression = evaluate_expression(nums, ops)
            if result == 24:
                return f"<<<{expression}={result}>>>"
    return "No solution found."

# Define the numbers
numbers = [3, 7, 7, 7]

# Find and print the expression that results in 24
print(find_expression_to_24(numbers))