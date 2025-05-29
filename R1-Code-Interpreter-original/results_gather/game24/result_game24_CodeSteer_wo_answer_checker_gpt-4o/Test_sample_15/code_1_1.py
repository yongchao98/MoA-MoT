from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Evaluate the expression with the given numbers and operations
    try:
        # Try different parenthesizations
        expr1 = f"(({nums[0]} {ops[0]} {nums[1]}) {ops[1]} ({nums[2]} {ops[2]} {nums[3]}))"
        expr2 = f"(({nums[0]} {ops[0]} ({nums[1]} {ops[1]} {nums[2]})) {ops[2]} {nums[3]})"
        expr3 = f"({nums[0]} {ops[0]} (({nums[1]} {ops[1]} {nums[2]}) {ops[2]} {nums[3]}))"
        expr4 = f"({nums[0]} {ops[0]} ({nums[1]} {ops[1]} ({nums[2]} {ops[2]} {nums[3]})))"
        expr5 = f"(({nums[0]} {ops[0]} {nums[1]}) {ops[1]} {nums[2]} {ops[2]} {nums[3]})"
        
        for expr in [expr1, expr2, expr3, expr4, expr5]:
            if eval(expr) == 24:
                return expr
    except ZeroDivisionError:
        pass
    return None

def find_expression_to_24(numbers):
    # Try all permutations of the numbers
    for nums in permutations(numbers):
        # Try all combinations of operations
        for ops in product("+-*/", repeat=3):
            result = evaluate_expression(nums, ops)
            if result:
                return result
    return "No solution found"

# Input numbers
numbers = [5, 5, 11, 12]

# Find and print the expression that evaluates to 24
solution = find_expression_to_24(numbers)
print(f"<<<{solution}>>>")