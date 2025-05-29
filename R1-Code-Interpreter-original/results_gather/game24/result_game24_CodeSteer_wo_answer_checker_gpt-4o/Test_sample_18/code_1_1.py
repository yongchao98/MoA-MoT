from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Evaluate the expression based on the current permutation of numbers and operations
    try:
        # Construct the expression as a string
        expression = f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]}))"
        # Evaluate the expression
        if eval(expression) == 24:
            return expression
    except ZeroDivisionError:
        return None
    return None

def find_expression_to_24(numbers):
    # Generate all permutations of the numbers
    for nums in permutations(numbers):
        # Generate all combinations of operations
        for ops in product("+-*/", repeat=3):
            # Evaluate the expression
            result = evaluate_expression(nums, ops)
            if result:
                return result
    return None

# Input numbers
numbers = [3, 7, 8, 8]

# Find the expression that evaluates to 24
expression = find_expression_to_24(numbers)

# Print the result in the required format
if expression:
    print(f"<<<{expression}>>>")
else:
    print("No solution found.")