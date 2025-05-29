from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Evaluate the expression based on the current permutation of numbers and operations
    try:
        result = nums[0]
        expression = str(nums[0])
        for i in range(3):
            if ops[i] == '+':
                result += nums[i + 1]
            elif ops[i] == '-':
                result -= nums[i + 1]
            elif ops[i] == '*':
                result *= nums[i + 1]
            elif ops[i] == '/':
                if nums[i + 1] == 0:
                    return None, None
                result /= nums[i + 1]
            expression += f" {ops[i]} {nums[i + 1]}"
        return result, expression
    except ZeroDivisionError:
        return None, None

def find_expression_to_24(numbers):
    # Generate all permutations of the numbers
    for nums in permutations(numbers):
        # Generate all combinations of operations
        for ops in product('+-*/', repeat=3):
            result, expression = evaluate_expression(nums, ops)
            if result == 24:
                return expression
    return None

# Input numbers
numbers = [3, 7, 10, 10]

# Find and print the expression that evaluates to 24
expression = find_expression_to_24(numbers)
if expression:
    print(f"<<<{expression}>>>")
else:
    print("No solution found.")