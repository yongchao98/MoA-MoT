from itertools import permutations
import operator

# Define the target value
TARGET = 24

# Define the operations
operations = {
    '+': operator.add,
    '-': operator.sub,
    '*': operator.mul,
    '/': operator.truediv
}

def evaluate_expression(nums, ops):
    """Evaluate the expression given numbers and operations."""
    try:
        result = nums[0]
        for i in range(3):
            result = ops[i](result, nums[i + 1])
        return result
    except ZeroDivisionError:
        return None

def find_expression(nums):
    """Find an expression that evaluates to the target value."""
    for num_perm in permutations(nums):
        for ops in permutations(operations.values(), 3):
            if evaluate_expression(num_perm, ops) == TARGET:
                # Construct the expression string
                expr = f"({num_perm[0]} {list(operations.keys())[list(operations.values()).index(ops[0])]} {num_perm[1]})"
                expr = f"({expr} {list(operations.keys())[list(operations.values()).index(ops[1])]} {num_perm[2]})"
                expr = f"({expr} {list(operations.keys())[list(operations.values()).index(ops[2])]} {num_perm[3]})"
                return expr
    return None

# Input numbers
numbers = [3, 3, 9, 9]

# Find and print the expression
expression = find_expression(numbers)
if expression:
    print(f"<<<{expression} = {TARGET}>>>")
else:
    print("No solution found.")