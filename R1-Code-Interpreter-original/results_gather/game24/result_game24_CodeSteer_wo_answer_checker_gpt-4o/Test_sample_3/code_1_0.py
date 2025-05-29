from itertools import permutations

def evaluate_expression(nums, ops):
    # Evaluate the expression based on the current permutation of numbers and operations
    result = nums[0]
    for i in range(3):
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

def find_expression(nums):
    operations = ['+', '-', '*', '/']
    for num_perm in permutations(nums):
        for ops in permutations(operations, 3):
            if evaluate_expression(num_perm, ops) == 24:
                return f"<<<({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}{num_perm[2]}{ops[2]}{num_perm[3]}>>>"
    return "No solution found"

# Input numbers
numbers = [2, 2, 2, 12]

# Find and print the expression
print(find_expression(numbers))