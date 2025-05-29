from itertools import permutations, product

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
                return None  # Avoid division by zero
            result /= nums[i + 1]
    return result

def find_expression(nums):
    # Try all permutations of numbers
    for num_perm in permutations(nums):
        # Try all combinations of operations
        for ops in product('+-*/', repeat=3):
            if evaluate_expression(num_perm, ops) == 24:
                # Format the expression as a string
                expression = f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})"
                return f"<<<{expression}>>>"
    return "No solution found"

# Input numbers
numbers = [1, 9, 12, 12]
print(find_expression(numbers))