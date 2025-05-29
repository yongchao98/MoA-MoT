from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Try different parenthesizations
    expressions = [
        f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}",
        f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}",
        f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})",
        f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]})",
        f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]}))"
    ]
    
    for expr in expressions:
        try:
            if abs(eval(expr) - 24) < 1e-6:  # Check if the result is 24
                return expr
        except ZeroDivisionError:
            continue
    return None

def find_expression_to_24(numbers):
    operators = ['+', '-', '*', '/']
    for nums in permutations(numbers):
        for ops in product(operators, repeat=3):
            result = evaluate_expression(nums, ops)
            if result:
                return f"<<<{result}>>>"
    return "No solution found"

numbers = [3, 6, 11, 11]
print(find_expression_to_24(numbers))