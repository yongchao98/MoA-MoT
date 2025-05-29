import itertools

def evaluate_expression(nums, ops):
    """Evaluate the expression with given numbers and operations."""
    try:
        # Try different parenthesis placements
        expressions = [
            f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}",
            f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}",
            f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})",
            f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]})",
            f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]}))"
        ]
        for expr in expressions:
            if eval(expr) == 24:
                return expr
    except ZeroDivisionError:
        pass
    return None

def find_expression_for_24(numbers):
    """Find an expression that evaluates to 24 using all numbers and basic operations."""
    operations = ['+', '-', '*', '/']
    for nums in itertools.permutations(numbers):
        for ops in itertools.product(operations, repeat=3):
            result = evaluate_expression(nums, ops)
            if result:
                return f"<<<{result}>>>"
    return "No solution found."

# Example usage
numbers = [6, 6, 12, 12]
print(find_expression_for_24(numbers))