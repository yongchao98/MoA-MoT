from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Evaluate the expression based on the current permutation of numbers and operations
    try:
        result = nums[0]
        for i in range(3):
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

def find_expression_to_24(numbers):
    # Generate all permutations of the numbers
    for num_perm in permutations(numbers):
        # Generate all combinations of operations
        for ops in product('+-*/', repeat=3):
            # Evaluate the expression
            if evaluate_expression(num_perm, ops) == 24:
                # Format the expression as a string
                expression = f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})"
                return f"<<<{expression}>>>"
    return "No solution found"

# Define the numbers
numbers = [1, 1, 2, 11]

# Find and print the expression that results in 24
print(find_expression_to_24(numbers))