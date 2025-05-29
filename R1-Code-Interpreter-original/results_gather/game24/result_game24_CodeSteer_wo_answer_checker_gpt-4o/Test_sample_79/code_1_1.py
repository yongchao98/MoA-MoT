from itertools import permutations
import operator

# Define the operations we can use
ops = {
    '+': operator.add,
    '-': operator.sub,
    '*': operator.mul,
    '/': operator.truediv
}

# Function to evaluate an expression
def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Recursive function to try all combinations of operations
def find_expression(nums, target=24):
    if len(nums) == 1:
        if abs(nums[0] - target) < 1e-6:  # Check if the result is close to 24
            return str(nums[0])
        return None

    for i in range(len(nums) - 1):
        for op in ops:
            left = nums[:i] + [f"({nums[i]}{op}{nums[i+1]})"] + nums[i+2:]
            result = evaluate_expression(left[0])
            if result is not None:
                expr = find_expression([result] + left[1:], target)
                if expr:
                    return expr
    return None

# Main function to find the solution
def solve_24(numbers):
    for perm in permutations(numbers):
        expr = find_expression(list(map(str, perm)))
        if expr:
            return f"<<<{expr}>>>"
    return "No solution found"

# Input numbers
numbers = [2, 6, 10, 12]

# Find and print the solution
print(solve_24(numbers))