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
def find_expression(nums, exprs, target=24):
    if len(nums) == 1:
        if abs(nums[0] - target) < 1e-6:  # Check if the result is close to 24
            return exprs[0]
        return None

    for i in range(len(nums) - 1):
        for op_symbol, op_func in ops.items():
            # Apply the operation to the first two numbers
            try:
                result = op_func(nums[i], nums[i+1])
            except ZeroDivisionError:
                continue

            # Construct the new expression
            new_expr = f"({exprs[i]}{op_symbol}{exprs[i+1]})"
            new_nums = nums[:i] + [result] + nums[i+2:]
            new_exprs = exprs[:i] + [new_expr] + exprs[i+2:]

            # Recursively find the expression
            expr = find_expression(new_nums, new_exprs, target)
            if expr:
                return expr
    return None

# Main function to find the solution
def solve_24(numbers):
    for perm in permutations(numbers):
        exprs = list(map(str, perm))
        expr = find_expression(list(perm), exprs)
        if expr:
            return f"<<<{expr}>>>"
    return "No solution found"

# Input numbers
numbers = [2, 6, 10, 12]

# Find and print the solution
print(solve_24(numbers))