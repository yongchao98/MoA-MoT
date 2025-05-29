from itertools import permutations
from sympy import symbols, simplify

def evaluate_expression(nums, ops):
    # Create symbolic variables for the numbers
    a, b, c, d = symbols('a b c d')
    # Create the expression using the given operations
    expr = f"(({a}{ops[0]}{b}){ops[1]}({c}{ops[2]}{d}))"
    # Substitute the numbers into the expression
    expr = expr.replace('a', str(nums[0])).replace('b', str(nums[1]))
    expr = expr.replace('c', str(nums[2])).replace('d', str(nums[3]))
    # Simplify the expression and check if it equals 24
    return simplify(expr) == 24, expr

def find_solution(numbers):
    # Define possible operations
    operations = ['+', '-', '*', '/']
    # Generate all permutations of the numbers
    for num_perm in permutations(numbers):
        # Generate all combinations of operations
        for ops in permutations(operations, 3):
            # Evaluate the expression
            is_valid, expr = evaluate_expression(num_perm, ops)
            if is_valid:
                return expr
    return None

# Given numbers
numbers = [5, 6, 11, 13]
solution = find_solution(numbers)

# Print the final answer in the required format
if solution:
    print(f"<<<{solution}>>>")
else:
    print("No solution found.")