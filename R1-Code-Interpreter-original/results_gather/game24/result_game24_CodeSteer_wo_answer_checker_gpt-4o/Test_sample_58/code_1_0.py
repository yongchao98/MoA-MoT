from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Evaluate the expression with the given numbers and operations
    try:
        result = eval(f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}")
        return result
    except ZeroDivisionError:
        return None

def find_expression_for_24(numbers):
    # Generate all permutations of the numbers
    for nums in permutations(numbers):
        # Generate all combinations of operations
        for ops in product("+-*/", repeat=3):
            # Evaluate the expression
            result = evaluate_expression(nums, ops)
            if result == 24:
                # Format the expression for output
                expression = f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}"
                return f"<<<{expression}>>>"
    return "No solution found."

# Test the function with the given input
numbers = [3, 7, 13, 13]
print(find_expression_for_24(numbers))