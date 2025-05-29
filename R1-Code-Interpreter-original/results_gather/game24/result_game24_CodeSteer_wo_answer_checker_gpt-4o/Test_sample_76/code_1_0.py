from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Define all possible expressions with different groupings
    expressions = [
        f"(({nums[0]} {ops[0]} {nums[1]}) {ops[1]} {nums[2]}) {ops[2]} {nums[3]}",
        f"({nums[0]} {ops[0]} ({nums[1]} {ops[1]} {nums[2]})) {ops[2]} {nums[3]}",
        f"{nums[0]} {ops[0]} (({nums[1]} {ops[1]} {nums[2]}) {ops[2]} {nums[3]})",
        f"{nums[0]} {ops[0]} ({nums[1]} {ops[1]} ({nums[2]} {ops[2]} {nums[3]}))",
        f"({nums[0]} {ops[0]} {nums[1]}) {ops[1]} ({nums[2]} {ops[2]} {nums[3]})"
    ]
    
    for expression in expressions:
        try:
            # Evaluate the expression
            result = eval(expression)
            # Check if the result is 24 and is an integer
            if result == 24 and result == int(result):
                return expression
        except ZeroDivisionError:
            # Ignore division by zero errors
            pass
    return None

def find_expression_to_24(numbers):
    # Define possible operations
    operations = ['+', '-', '*', '/']
    # Generate all permutations of the numbers
    for nums in permutations(numbers):
        # Generate all combinations of operations
        for ops in product(operations, repeat=3):
            # Evaluate the expression
            expression = evaluate_expression(nums, ops)
            if expression:
                return expression
    return "No solution found"

numbers = [6, 8, 11, 11]
solution = find_expression_to_24(numbers)
print(f"<<<{solution}>>>")