from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Construct the expression as a string
    expr = f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})"
    try:
        # Evaluate the expression
        if eval(expr) == 24:
            return expr
    except ZeroDivisionError:
        pass
    return None

def find_expression_to_24(numbers):
    # Generate all permutations of the numbers
    for nums in permutations(numbers):
        # Generate all combinations of operations
        for ops in product("+-*/", repeat=3):
            # Evaluate the expression
            expr = evaluate_expression(nums, ops)
            if expr:
                return expr
    return "No solution found"

# Define the numbers
numbers = [1, 2, 11, 11]

# Find and print the expression that evaluates to 24
solution = find_expression_to_24(numbers)
print(f"<<<{solution}>>>")